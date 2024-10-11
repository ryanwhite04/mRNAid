import json
from flask import Flask, request, render_template, jsonify
from flask_socketio import SocketIO, emit, join_room, leave_room
from flask_cors import CORS
import pickle
from tasks import optimization_evaluation_task, arwa_generator_task
from common.utils.Exceptions import EmptySequenceError, SequenceLengthError, NoGCError, EntropyWindowError, \
    NumberOfSequencesError, WrongCharSequenceError, RangeError, SpeciesError
from common.utils.Logger import MyLogger
from common.utils.RequestParser import RequestParser
from common.arw_mrna.src import protein, awalk, vienna, objective_functions as objectives
import python_codon_tables as pct
from os import getcwd, path, environ
from notify import send_email
app = Flask(__name__)
app.config['SECRET_KEY'] = environ.get('SECRET_KEY', 'default-secret-key')
socketio = SocketIO(app)
CORS(app)

# Setting up a logger
logger = MyLogger(__name__)


@app.route('/codon_table', methods=['GET'])
def codon_table():
    return render_template('codon_table.html')

@app.route('/codon_table', methods=['POST'])
def codon_table_post():
    # receive form input, get taxid 
    taxid = request.form['taxid']
    cd = pct.get_codons_table(taxid)
    if cd is None:
        return jsonify({'error': 'Invalid taxid'}), 400
    # convert the cd dictionary to a json object and return it
    return jsonify(cd)
    
@app.route('/test_email', methods=['GET'])
def test_email_form():
    return render_template('test_email.html')

@app.route('/test_email', methods=['POST'])
def test_email():
    try:
        # Extract data from the form
        subject = request.form["subject"]
        body = request.form["message"]  # Note: You had "body" before, but in the form, it's named "message"
        recipient = request.form["email"]  # Note: You had "recipient" before, but in the form, it's named "email"
        
        # Send the email using your notify.py function
        send_email(subject, body, recipient)
        
        return jsonify({'message': 'Email sent successfully!'}), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/custom_codon_table/<taxID>', methods=['GET'])
def custom_codon_table(taxID):
    try:
        print(taxID)
            # LOAD ONE TABLE BY TAXID (it will get it from the internet if it is not
        # in the builtin tables)
        table = pct.get_codons_table(taxID)
        return jsonify(table), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500
    
@socketio.on('task_progress')
def handle_task_progress(data):
    room = data['task_id']
    step = data['result']
    socketio.emit('task_update', step, room=room)

@app.route('/task/<task_id>', methods=['GET'])
def get_task(task_id):
    task = arwa_generator_task.AsyncResult(task_id)
    data = None
    if task.status == "SUCCESS":
        data = task.get()
    return render_template('task.html', task=task, data=data)

@app.route('/arwa_sync', methods=['GET'])
def arwa_sync_form():
    return render_template('arwa_sync.html')

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
        elif path.exists('common/arw_mrna/codon_tables'):
            base_dir = 'common/arw_mrna/codon_tables'
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
        return jsonify({
            'cds': result.cds,
            'fitness': result.fitness,
            'cai': cai,
            'aup': aup,
            'efe': efe
        }), 200
    except Exception as e:
        logger.error(f'Error handling ARWA request: {str(e)}', exc_info=True)
        return jsonify({'error': str(e)}), 500

@app.route('/arwa_websocket', methods=['GET'])
def arwa_websocket():
    return render_template('arwa_websocket.html')

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

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
@socketio.on('join')
def on_join(data):
    room = data['room']
    join_room(room)

@socketio.on('leave')
def on_leave(data):
    room = data['room']
    leave_room(room)

@socketio.on('arwa_websocket_celery')
def handle_arwa_websocket_celery(args):
    def on_message(body):
        status = body["status"]
        if status in ["final", "initial", "progress", "new_cds"]:
            emit('arwa_sync_progress', body["result"])
    try:
        print('Socket IO args')
        print(args)
        args["verbose"] = True
        task = arwa_generator_task.delay(args)
        join_room(task.id)
        # task.get(on_message=on_message, propagate=False)
        
    except Exception as e:
        logger.error(f'Error handling ARWA request: {str(e)}', exc_info=True)
        emit('arwa_sync_error', {'error': str(e)})
        
@socketio.on('arwa_websocket')
def handle_arwa_sync(args):
    try:
        cai_threshold = float(args["cai_threshold"])
        cai_exp_scale = float(args["cai_exp_scale"])
        stability = args["stability"]
        verbose = True
        freq_table_path = args["freq_table_path"]
        taxid = freq_table_path if freq_table_path.endswith('.txt') else int(freq_table_path)
        freq_table = protein.CodonFrequencyTable(taxid)
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
                        # only send progress updates if they are a multiple of 10
            if update["type"] == "progress" and update["step"] % 10 != 0:
                continue
            update = { **update, "stability": stability, "cai_threshold": cai_threshold, "cai_exp_scale": cai_exp_scale}
            emit('arwa_sync_progress', update)
            # Send an email when the task is complete
            if update["type"] == "final" and args["email"]:
                print("Sending email")
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
        task = arwa_generator_task.delay(parameters)
        # TODO, this is blocking.
        # print(task.get(on_message=print, propagate=False))
        return jsonify({'task_id': task.id}), 200
    except Exception as e:
        logger.error(f'Error handling ARWA request: {str(e)}', exc_info=True)
        return jsonify({'error': str(e)}), 500

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
        return jsonify({'task_id': task.id})


@app.route('/api/v1/status/<task_id>', methods=['GET'])
def status(task_id: str) -> str:
    """
    Get task status based on task ID
    :param task_id: string
    :return: JSON with task ID, task state, message and task result, if any.
    """
    task = optimization_evaluation_task.AsyncResult(task_id)
    if task.state == 'PENDING':
        return jsonify(
            {"message": "Task is in progress... Lay back or get some coffee, I'm working hard on your request!",
             "task_id": task.id,
             "state": task.state,
             "data": None}), 202
    elif task.state == 'FAILURE':
        logger.info(10 * '#' + 'END OF PROCESSING THE REQUEST: FAILURE' + 10 * '#')
        return jsonify(
            {
                "message": f"Ooops, something went wrong! Please contact administrator and provide him/her your task id: {task.id}",
                "task_id": task.id,
                "state": task.state,
                "data": None}), 500
    elif task.state == 'SUCCESS':
        logger.info(10 * '#' + 'END OF PROCESSING THE REQUEST: SUCCESS' + 10 * '#')
        return jsonify(
            {"message": "Done!",
             "task_id": task.id,
             "state": task.state,
             "data": json.loads(task.get())})

application = app