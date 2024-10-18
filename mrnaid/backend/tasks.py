import json
import numpy as np
from billiard import Pool
from celery import Celery
from celery.signals import task_prerun, task_postrun
from socketio import Client
from common.Evaluation import Evaluation
from common.OptimizationProblems import initialize_optimization_problem
from common.OptimizationTask import optimization_task
from common.utils.Datatypes import OptimizationParameters
from common.utils.Logger import MyLogger
from common.arw_mrna.src import protein, awalk, objective_functions as objectives
from notify import send_email
from datetime import datetime, timedelta
from config import PRIVATE_URL, PUBLIC_URL, NUMBER_OF_ATTEMPTS

logger = MyLogger(__name__)
celery = Celery('tasks')
celery.config_from_object('celery_config')
sio = Client()

def init_pool():
    np.random.seed()

@task_prerun.connect
def task_prerun_handler(sender=None, **kwargs):
    logger.info(F"{PRIVATE_URL=}")
    if not sio.connected:
        sio.connect(PRIVATE_URL)

@task_postrun.connect
def task_postrun_handler(sender=None, **kwargs):
    if sio.connected:
        sio.disconnect()

@celery.task(bind=True)
def arwa_task(self, args: dict) -> str:
    logger.info(10 * '#' + 'STARTING PROCESSING THE REQUEST' + 10 * '#')
    logger.info(f"Received request with args: {args}")
    task_id = self.request.id
    sio.call('join', {'room': task_id})
    cai_threshold = args["cai_threshold"]
    cai_exp_scale = args["cai_exp_scale"]
    # Construct the absolute path to the freq_table file
    stability = args["stability"]
    verbose = args["verbose"]
    freq_table_path = args["freq_table_path"]
    seed = args.get("seed")
    taxid = freq_table_path if freq_table_path.endswith('.txt') else int(freq_table_path)
    freq_table = protein.CodonFrequencyTable(taxid)
    # Get obj function
    obj = {
        'aup': objectives.make_cai_and_aup_obj,
        'efe': objectives.make_cai_and_efe_obj,
        'none': objectives.make_cai_threshold_obj
    }[stability](objectives.CAIThresholdObjectiveConfig(
        freq_table,
        cai_threshold,
        cai_exp_scale,
        verbose=verbose
    ))
    config = awalk.WalkConfig(
        args["aa_seq"],
        freq_table,
        obj,
        args["steps"],
        seed=seed,
    )
    # Create walk config]
    result: awalk.WalkResult
    for result in awalk.adaptive_random_walk(config):
        sio.call('task_progress', {'task_id': task_id, 'result': result.to_dict()}, timeout=2)
    if args["email"]:
        expiry = datetime.now() + timedelta(days=30)
        subject = "Task Complete"
        body = f"""
            View the result at {PUBLIC_URL}/task/{task_id}
            This link will expire at {expiry.strftime("%Y-%m-%d %H:%M:%S")}
        """
        send_email(subject, body, args["email"])
        logger.info(f"Email sent to {args['email']}")
    # wait for the last message to be sent
    
    sio.call('leave', {'room': task_id})
    logger.info(f"Final result: {result.to_dict()}")
    return json.dumps(result.to_dict())

@celery.task()
def optimization_evaluation_task(parameters: dict) -> str:
    """
    Optimize sequences using Optimization class, evaluate results using Evaluation class
    :param parameters: Dictionary created from OptimizationParameters class
    :return: JSON with the response dictionary containing both input and optimized sequences and their evaluation
    """

    parameters = OptimizationParameters.parse_obj(parameters)
    optimization_problems = [initialize_optimization_problem(parameters) for _ in range(parameters.number_of_sequences)]

    pool = Pool(initializer=init_pool)
    result = pool.map(optimization_task, optimization_problems)

    if len(set(result)) < parameters.number_of_sequences:
        logger.info('Less sequences than required!')
        attempt_count = 0
        while len(set(result)) < parameters.number_of_sequences and attempt_count < NUMBER_OF_ATTEMPTS:
            logger.info('Performing additional optimization attempt')
            result.append(optimization_task(initialize_optimization_problem(parameters)))
            attempt_count += 1

            if len(set(result)) == parameters.number_of_sequences:
                logger.info(f'Correct number achieved! Stopping...')
                break
            elif (attempt_count == parameters.number_of_sequences) and (
                    len(set(result)) < parameters.number_of_sequences):
                logger.info('Not able to provide required number of results. Stopping...')
                break
            else:
                pass

    evaluator = Evaluation(list(set(result)), parameters)
    final_response = evaluator.get_evaluation()
    final_response["input_parameters"] = parameters.dict()
    return json.dumps(final_response)
