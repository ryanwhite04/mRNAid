# following commands need to work from backend folder

make gunicorn
make redis
make worker
make flower
make flask
make test

gunicorn app
celery worker -A tasks --loglevel=info
celery flower -A tasks --port=5555
flask --app app run --host=0.0.0.0 # wont' work because of the way flask is imported
python -m app # this will work instead because of the way flask is imported
# flask imports app from app.py, so we can run it as a module
python -m unittest discover
redis-server, can use ./redis/src/red


gunicorn app:app
celery worker -A tasks --loglevel=info
celery flower -A tasks --port=5555
flask --app app run --host=0.0.0.0
python -m unittest discover
redis-server, can use ./redis/src/red