import json
import os
import numpy as np
from billiard import Pool
from celery import Celery
from celery.signals import task_prerun, task_postrun
from socketio import Client
import pickle
from common.Evaluation import Evaluation
from common.OptimizationProblems import initialize_optimization_problem
from common.OptimizationTask import optimization_task
from common.utils.Datatypes import OptimizationParameters
from common.utils.Logger import MyLogger
from common.arw_mrna.src import protein, awalk, vienna, objective_functions as objectives
from notify import send_email
from datetime import datetime, timedelta

# Setting up logger
logger = MyLogger(__name__)

NUMBER_OF_ATTEMPTS = 3
PRIVATE_HOST = os.getenv('PRIVATE_HOST', "http://flask:5000")
PUBLIC_HOST = os.getenv('PUBLIC_HOST', "http://localhost:5000")
celery = Celery('tasks')
celery.config_from_object('celery_config')
sio = Client()

# @sio.event
# def connect():
#     print("Celery worker connected to WebSocket")

# @sio.event
# def disconnect():
#     print("Celery worker disconnected from WebSocket")

@task_prerun.connect
def task_prerun_handler(sender=None, **kwargs):
    if not sio.connected:
        sio.connect(PRIVATE_HOST)

@task_postrun.connect
def task_postrun_handler(sender=None, **kwargs):
    if sio.connected:
        sio.disconnect()

def init_pool():
    np.random.seed()

# @celery.task
# def ARWA(seconds: int):
#     try:
#         # Log start of task
#         logger.info(f"Starting ARWA task with seconds: {seconds}")
        
#         # Simulate a long-running task
#         time.sleep(seconds)  # Task takes 'seconds' seconds to complete
        
#         # Log completion
#         logger.info("ARWA task completed successfully.")
#         response = {}
#         response["message"] = f"Task Completed after {seconds} seconds"
#         return json.dumps(response)
#     except Exception as e:
#         logger.error(f"Error in ARWA task: {str(e)}", exc_info=True)
#         raise e

@celery.task()
def ARWA(args: dict) -> str:
    print("hi")
    logger.info(10 * '#' + 'STARTING PROCESSING THE REQUEST' + 10 * '#')
    logger.info(f"Received request with args: {args}")
    logger.info("STARTING PROCESSING THE REQUEST")
    logger.info(f"Received request with args: {args}")
    try:
        cai_threshold = args["cai_threshold"]
        cai_exp_scale = args["cai_exp_scale"]
        # Construct the absolute path to the freq_table file
        base_path = os.path.dirname(os.path.abspath(__file__))
        freq_table_path = os.path.join(base_path, 'common', 'arw_mrna', 'codon_tables', args['freq_table_path'])
        logger.info(f"Looking for freq_table at {freq_table_path}")

        if not os.path.exists(freq_table_path):
            raise FileNotFoundError(f"The file does not exist: {freq_table_path}")

        # freq_table_path = f"mrnaid/backend/common/arw_mrna/codon_tables/{args['freq_table_path']}"
        stability = args["stability"]
        verbose = args["verbose"]
        freq_table = protein.CodonFrequencyTable(freq_table_path)
        obj_config = objectives.CAIThresholdObjectiveConfig(
            freq_table,
            cai_threshold,
            cai_exp_scale,
            verbose=verbose
        )
        # Get obj function
        if stability == 'aup':
            obj = objectives.make_cai_and_aup_obj(obj_config)
        elif stability == 'efe':
            obj = objectives.make_cai_and_efe_obj(obj_config)
        elif stability == 'none':
            obj = objectives.make_cai_threshold_obj(obj_config)
        init_cds = None
        load_path = args.get("load_path")
        if load_path is not None:
            with open(load_path, "rb") as f:
                init_cds = pickle.load(f).cds
        # Create walk config]
        aa_seq = args["aa_seq"]
        steps = args["steps"]
        walk_config = awalk.WalkConfig(
            aa_seq, freq_table, obj, steps, init_cds=init_cds, verbose=verbose)
        result = awalk.adaptive_random_walk(walk_config)
        cai = freq_table.codon_adaptation_index(result.cds)
        aup, efe = vienna.aup_and_efe(vienna.cds_to_rna(result.cds))
        logger.info(10 * '#' + 'END OF PROCESSING THE REQUEST: SUCCESS' + 10 * '#')
        logger.info("END OF PROCESSING THE REQUEST: SUCCESS")
        return json.dumps({
            'cds': result.cds,
            'fitness': result.fitness,
            'cai': cai,
            'aup': aup,
            'efe': efe
        })
    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=True)
        raise

@celery.task(bind=True)
def arwa_generator_task(self, args: dict) -> str:
    logger.info(10 * '#' + 'STARTING PROCESSING THE REQUEST' + 10 * '#')
    logger.info(f"Received request with args: {args}")
    task_id = self.request.id
    if not sio.connected: sio.connect(PRIVATE_HOST)
    sio.emit('join', {'room': task_id})
    cai_threshold = args["cai_threshold"]
    cai_exp_scale = args["cai_exp_scale"]
    # Construct the absolute path to the freq_table file
    stability = args["stability"]
    verbose = args["verbose"]
    freq_table_path = args["freq_table_path"]
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
    init_cds = None
    load_path = args.get("load_path")
    if load_path is not None:
        with open(load_path, "rb") as f:
            init_cds = pickle.load(f).cds
            
    # Create walk config]
    generator = awalk.adaptive_random_walk_generator(awalk.WalkConfig(
        args["aa_seq"],
        freq_table,
        obj,
        args["steps"],
        init_cds=init_cds,
        verbose=verbose
    ))
    for step in generator:
        # only send progress updates if they are a multiple of 10
        if step["type"] == "progress" and step["step"] % 10 != 0:
            continue
        result = {**step, "stability": stability, "cai_threshold": cai_threshold, "cai_exp_scale": cai_exp_scale}
        # self.update_state(state=step["type"], meta=step)
        sio.emit('task_progress', {'task_id': task_id, 'result': result})
    # sio.emit('task_progress', {'task_id': task_id, 'status': 'SUCCESS', 'result': step})

    if args["email"]:
        print("Sending email")
        expiry = datetime.now() + timedelta(days=30)
        subject = "Task Complete"
        body = f"""
            View the result at {PUBLIC_HOST}/task/{task_id}
            This link will expire at {expiry.strftime("%Y-%m-%d %H:%M:%S")}
        """
        send_email(subject, body, args["email"])
    sio.emit('leave', {'room': task_id})
    if sio.connected: sio.disconnect()
    return json.dumps(result)

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
