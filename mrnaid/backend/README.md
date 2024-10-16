
# Development

This project contains a submodule so remember to

```
git submodule update --init --recursive
```

# Basic Information

This app has 4 main components

## Flask
    This is the main server that serves the frontend and the api
    It also handles the celery tasks

    This can be replaced with *gunicorn* for production

## Worker (Celery Worker)

This is the task queue that runs the arwa algorithm
It is used to run the arwa algorithm in the background
It can be scaled to run on multiple servers

It uses redis as the message broker
It also has a flower dashboard for monitoring tasks

To run locally
```sh
celery -A tasks worker --loglevel=info
```

You can change the amount of workers with the -c flag
```sh
celery -A tasks worker -c 4 --loglevel=info
```

To run in docker
```sh
docker compose up worker
```

If you are running in docker and want more workers
You have to update the docker-compose.yml file
For now it's just default because that seems to work fine

## Redis
    
To run locally
```sh
redis-server
```
This requires installing redis locally obviously
There is some code in here for automatically pulling redis and running it from ./redis folder but it's outdated and didn't work for everyone, easier to just install normally or use docker

To run in docker
```sh
docker compose up redis
```

This is the message broker for celery
It is used to store the tasks and results

## Flower
    
To run locally
```sh
celery -A tasks flower
```

To run in docker
```sh
docker compose up flower
```
This is a dashboard for monitoring the tasks
It is useful for debugging tasks and monitoring the workers

It is available on port 5555 by default

You'll need to login with the FLOWER_USERNAME and FLOWER_PASSWORD from the .env file

# Directory Structure

All commands should be run from within this directory, mrnaid/backend
The other directories are for other features in the original codebase that we aren't using

- common/
    This is common utilities for the project
    it include arw_mrna, which has the arwa code
- docs/
    This is for original codebase, ignore this
- logs/
    Logs from workers are stored here in logs.log
    You can change this with the LOG_FILE environment variable
    Code for logger is in common/utils/Logger.py
- static/
    Static content for frontend, just index.js and index.css
    This includes styling for frontend and the code for the charts and websocket
- templates/
    Page templates in jinja2
    Only pages to worry about are index.html and codon_table.html
    These are all no similarly named routes in route.py
- tests/
    Tests for the project
    Currently only 2 tests that are automated
    Testing arwa algorithm and that the codon_table class is correct
- app.py
    Main file for the project
    This is where the flask app is created
    used for flask run
- celery_config.py
    Configuration for celery
    This includes redis connection urls and expiration date on tasks
    expiration is currently 30 days
- docker-compose.yml
    Docker compose file for the project
    This includes the redis, worker, flower, gunicorn, and flask services
    Don't use the docker-compose in the root directory
- Dockerfile
- environment.yml
    Environment file for conda
    ```conda env create -f environment.yml```
- notify.py
    This is the file that sends emails
    It uses the sendgrid api
    Can probably move this into common/utils later
- requirements.txt
    This is used by docker so that it doesn't require conda
    can create with
    ```conda list --format=freeze > requirements.txt```
- routes.py
    This is where all the routes are defined
    This includes the routes for the frontend and the api
- tasks.py
    This is where all the celery tasks are defined
    This includes the tasks for the arwa algorithm and the original mrnaid algorithm

# Environment Setup

.env should include the following

- SENDGRID_API_KEY:
  Can get this when you sign up on sendgrid
- SENDGRID_EMAIL_USERNAME
  get this when signing up on sendgrid
- REDIS_ADDRESS
    ip address of the redis server, if running locally this is "flask"
    if you running a worker on a separate server, you the ip address of the redis server in the private subnet
- REDIS_PASSWORD
    password of the redis server
    if you don't set this there will be a warning
- PRIVATE_HOST
    url of public flask server in subnet
    for instance, http://flask:5000
- PUBLIC_HOST
    url of public flask server
    for instance http://localhost:5000
- FLOWER_USERNAME
- FLOWER_PASSWORD
    Use these to login to the flower dashboard
- SECRET_KEY
    secret key for the flask app
    can be generated with
    ```python
    import os
    os.urandom(24)
    ```    

# Docker

To build the docker image, run the following command in the root directory of the project:
```
docker compose up redis
docker compose up flower
docker compose up worker
docker compose up flask
```

# Tests

Currently only 2 tests that are automated

Testing arwa algorithm and that the codon_table class is correct
As these 2 have changed a lot from the original arwa codebase
```
python -m unittest tests.test_arwa
python -m unittest tests.test_codon_table
```

the other tests are from previous codebase and weren't working when pulled
Haven't been fixed yet, you can delete them if you want, but they might be useful later on