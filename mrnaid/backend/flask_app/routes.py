import json

from flask import Flask, request, render_template
from tasks import optimization_evaluation_task, ARWA
from flask_socketio import SocketIO, emit
from flask_cors import CORS
from utils.Exceptions import EmptySequenceError, SequenceLengthError, NoGCError, EntropyWindowError, \
    NumberOfSequencesError, WrongCharSequenceError, RangeError, SpeciesError
from utils.Logger import MyLogger
from utils.RequestParser import RequestParser
import pickle
import protein
import awalk
import vienna
import objective_functions as objectives
from os import getcwd, path
import smtplib
from email.mime.text import MIMEText

app = Flask(__name__)
socketio = SocketIO(app)
CORS(app)

# Setting up a logger
logger = MyLogger(__name__)

# Determine if we are running from the project root or from within the flask app
cwd = getcwd()
if path.basename(cwd) == 'flask_app':
    project_root = path.abspath(path.join(cwd, '../../..'))
else:
    project_root = path.abspath(cwd)

def send_email(subject, body, recipient):
    
    # Use the email configuration from email_config.py
    
    from email_config import email_server, email_port, gmail_username, gmail_password
    sender = gmail_username
    msg = MIMEText(body)
    msg["From"] = sender
    msg["To"] = recipient
    msg["Subject"] = subject

    with smtplib.SMTP(email_server, email_port) as server:
        server.starttls()
        server.login(gmail_username, gmail_password)
        server.sendmail(sender, recipient, msg.as_string())

    print("Email sent successfully!")

@app.route('/arwa_sync', methods=['GET'])
def arwa_sync_form():
    return render_template('arwa_sync.html')

#     ap.add_argument('--steps', type=int, default=1000,
#                     help='Number of steps in the adaptive walk')
#     ap.add_argument('--aa_seq', type=str, required=True,
#                     help='Amino acid sequence to find a CDS for. A string using the standard one-letter code.')
#     ap.add_argument('--stability', type=str, choices=['aup', 'efe', 'none'], default='aup',
#                     help='Stability objective to use. Set to none to only use CAI')
#     ap.add_argument('--freq_table_path', type=str,
#                     default='../codon_tables/homosapiens.txt', help='Path to codon frequency table')
#     ap.add_argument('--cai_threshold', type=float, default=0.8,
#                     help='Objective function forces CAI to be at least this')
#     ap.add_argument('--cai_exp_scale', type=float, default=1.0,
#                     help='Scaling factor for CAI. Increase to make CAI more important')
#     ap.add_argument("--save_path", type=str, default=None,
#                     help="The path to save the result. Saved in pickle format.")
#     ap.add_argument("--verbose", action="store_true", help="Log all progress")
#     ap.add_argument("--load_path", type=str, default=None,
#                     help="Loads the initial CDS from the given pickle file. If not specified, the initial CDS is randomly generated.")
#     args = ap.parse_args()

#     # Load frequency table
#     freq_table = protein.CodonFrequencyTable(args.freq_table_path)

#     obj_config = objectives.CAIThresholdObjectiveConfig(
#         freq_table, args.cai_threshold, args.cai_exp_scale, verbose=args.verbose)

#     # Get obj function
#     if args.stability == 'aup':
#         obj = objectives.make_cai_and_aup_obj(obj_config)
#     elif args.stability == 'efe':
#         obj = objectives.make_cai_and_efe_obj(obj_config)
#     elif args.stability == 'none':
#         obj = objectives.make_cai_threshold_obj(obj_config)

#     # Load initial CDS if specified
#     init_cds = None
#     if args.load_path is not None:
#         # Read previous result as pickle
#         with open(args.load_path, "rb") as f:
#             init_cds = pickle.load(f).cds

#     # Create walk config
#     walk_config = awalk.WalkConfig(
#         args.aa_seq, freq_table, obj, args.steps, init_cds=init_cds, verbose=args.verbose)
@app.route('/api/v1/arwa_sync', methods=['POST'])
def awra_sync():
    logger.info(10 * '#' + 'NEW REQUEST' + 10 * '#' + 'SYNC')
    try:
        if request.content_type == 'application/json':
            args = request.get_json()
        else:
            args = request.form.to_dict()
            args["cai_threshold"] = float(args["cai_threshold"])
            args["cai_exp_scale"] = float(args["cai_exp_scale"])
            args["verbose"] = args["verbose"] == "true"
            args["steps"] = int(args["steps"])
        cai_threshold = args["cai_threshold"]
        cai_exp_scale = args["cai_exp_scale"]
        print("cwd: ", getcwd())
        # Determine the base directory for the codon tables
        if path.exists('mrnaid/backend/common/arw_mrna/codon_tables'):
            base_dir = 'mrnaid/backend/common/arw_mrna/codon_tables'
        elif path.exists('../common/arw_mrna/codon_tables'):
            base_dir = '../common/arw_mrna/codon_tables'
        else:
            raise FileNotFoundError("Cannot find the codon tables directory")

        # Construct the path to the frequency table
        freq_table_path = path.join(base_dir, args['freq_table_path'])
        print(freq_table_path)
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
        return json.dumps({
            'cds': result.cds,
            'fitness': result.fitness,
            'cai': cai,
            'aup': aup,
            'efe': efe
        }), 200
    except Exception as e:
        logger.error(f'Error handling ARWA request: {str(e)}', exc_info=True)
        return json.dumps({'error': str(e)}), 500

@app.route('/arwa_websocket', methods=['GET'])
def arwa_websocket():
    return render_template('arwa_websocket.html')

def format_args(args):
    # for printing args nicely in email
    return f"""
        CAI Threshold: {args["cai_threshold"]}
        CAI Exp Scale: {args["cai_exp_scale"]}
        Species: {args["freq_table_path"]}
        Stability: {args["stability"]}
        Amino Acids Sequence: {args["aa_seq"]}
        Steps: {args["steps"]}
        
    """
@socketio.on('arwa_websocket')
def handle_arwa_sync(args):
    try:
        cai_threshold = float(args["cai_threshold"])
        cai_exp_scale = float(args["cai_exp_scale"])
        stability = args["stability"]
        verbose = True
        print("cwd: ", getcwd())
        # Determine the base directory for the codon tables
        if path.exists('mrnaid/backend/common/arw_mrna/codon_tables'):
            base_dir = 'mrnaid/backend/common/arw_mrna/codon_tables'
        elif path.exists('../common/arw_mrna/codon_tables'):
            base_dir = '../common/arw_mrna/codon_tables'
        else:
            raise FileNotFoundError("Cannot find the codon tables directory")

        # Construct the path to the frequency table
        freq_table_path = path.join(base_dir, args['freq_table_path'])
        print(freq_table_path)
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
        if load_path:
            with open(load_path, "rb") as f:
                init_cds = pickle.load(f).cds
            
        # Create walk config
        aa_seq = args["aa_seq"]
        steps = int(args["steps"])
        walk_config = awalk.WalkConfig(
            aa_seq, freq_table, obj, steps, init_cds=init_cds, verbose=verbose)

        for update in awalk.adaptive_random_walk_generator(walk_config):
            emit('arwa_sync_progress', update)
            # Send an email when the task is complete
            if update["type"] == "final" and args["email"]:
                subject = "Task Complete"
                body = f"{format_args(args)}\n\nCDS: {update['cds']}\nFitness: {update['fitness']}"
                send_email(subject, body, args["email"])
                
    except Exception as e:
        logger.error(f'Error handling ARWA request: {str(e)}', exc_info=True)
        emit('arwa_sync_error', {'error': str(e)})
    
@app.route('/api/v1/arwa', methods=['POST'])
def arwa():
    logger.info(10 * '#' + 'NEW REQUEST' + 10 * '#')
    try:
        parameters = json.loads(request.get_data())
        logger.debug('Sequences are received from Parser')
        task = ARWA.delay(10)
        return json.dumps({'task_id': task.id}), 200
    except Exception as e:
        logger.error(f'Error handling ARWA request: {str(e)}', exc_info=True)
        return json.dumps({'error': str(e)}), 500

# Defining the routes
@app.route('/api/v1/optimize', methods=['POST'])
def optimization():
    """
    Parse optimization task parameters from the request, create a celery task for optimization.
    :return: JSON with task ID
    """
    # Processing the input
    logger.info(10 * '#' + 'NEW REQUEST' + 10 * '#')
    parser = RequestParser(request)

    # -- Getting input sequences
    try:
        parameters = parser.parse()
    except (EmptySequenceError, SequenceLengthError, NoGCError, EntropyWindowError, NumberOfSequencesError,
            WrongCharSequenceError,
            RangeError, SpeciesError) as e:
        logger.error(e.message)
        return e.message, 500
    except KeyError as e:
        logger.error('Error with json structure for the input sequences', exc_info=True)
        return 'Problem with json structure for the input sequences. Contact administrator. {}'.format(str(e)), 500
    else:
        logger.debug('Sequences are received from Parser')

    # Task optimization/evaluation:
    try:
        task = optimization_evaluation_task.apply_async(args=[parameters.dict()],
                                                        kwargs={})
    except Exception as e:
        import traceback
        traceback.print_exc()
        logger.error(f'Something went wrong when sending task to the queue: {e}')
        return str(e), 500
    else:
        return json.dumps({'task_id': task.id})


@app.route('/api/v1/status/<task_id>', methods=['GET'])
def status(task_id: str) -> str:
    """
    Get task status based on task ID
    :param task_id: string
    :return: JSON with task ID, task state, message and task result, if any.
    """
    task = optimization_evaluation_task.AsyncResult(task_id)
    if task.state == 'PENDING':
        return json.dumps(
            {"message": "Task is in progress... Lay back or get some coffee, I'm working hard on your request!",
             "task_id": task.id,
             "state": task.state,
             "data": None}), 202
    elif task.state == 'FAILURE':
        logger.info(10 * '#' + 'END OF PROCESSING THE REQUEST: FAILURE' + 10 * '#')
        return json.dumps(
            {
                "message": f"Ooops, something went wrong! Please contact administrator and provide him/her your task id: {task.id}",
                "task_id": task.id,
                "state": task.state,
                "data": None}), 500
    elif task.state == 'SUCCESS':
        logger.info(10 * '#' + 'END OF PROCESSING THE REQUEST: SUCCESS' + 10 * '#')
        return json.dumps(
            {"message": "Done!",
             "task_id": task.id,
             "state": task.state,
             "data": json.loads(task.get())})
