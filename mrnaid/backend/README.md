
# Development

This project contains a submodule so remember to

```
git submodule update --init --recursive
```

# Basic Information

This app has 5 services
Before starting gunicorn, flask, or worker, make sure you have the conda environment activated if using conda

You also need redis running for the worker, flask, or flower to run, and the worker needs flask running

## Flask

This is the main server that serves the frontend and the api
It also handles the celery tasks

This can be replaced with *gunicorn* for production

To run locally
```sh
flask run
```

To run in docker
```sh
docker compose up flask
```

## Gunicorn

This is a production server that can be used to run the flask app
If you use this instead remember to update the "private_host" and "public_host" in the .env file
So that they point to the gunicorn port instead of the flask port

To run locally
```sh
gunicorn -w 4 -b
```

To run in docker
```sh
docker compose up gunicorn
```

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

By default, any local .env folder will be used to set environmet variables
So make sure these also work for local development if not running in docker
And in particular, PRIVATE_HOST is not flask

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
    ```python -m pip list --format=freeze > requirements.txt```
- routes.py
    This is where all the routes are defined
    This includes the routes for the frontend and the api
- tasks.py
    This is where all the celery tasks are defined
    This includes the tasks for the arwa algorithm and the original mrnaid algorithm

# Environment Setup

In order to run the docker containers you need a .env file in the mrnaid/backend folder with the following variables:

```bash
# Sendgrid api details
# get these when sign up on sendgrid
# If not set, app will work but emails won't be sent
SENDGRID_API_KEY=<api key>
SENDGRID_EMAIL_USERNAME=<username>

# Redis details
# ip address of the redis server, if running locally this is "redis"
# if you running a worker on a separate server, you the ip address of the redis server in the private subnet
REDIS_ADDRESS=<address>

# password of the redis server
# if you don't set this there will be a warning, and you can't bind to 0.0.0.0
REDIS_PASSWORD=<password>

# Flower details
# Use these to log into flower
# using credentials because it will be exposed on a digital ocean droplet
FLOWER_USERNAME=<username>
FLOWER_PASSWORD=<password>

# Private server details
# These are the server details the worker will use for it's websocket connection
# If running locally will just be localhost 5000
# If running in docker will be flask 5000 or gunicorn 5000
# If running on digital ocean, use the subnet ip address of the droplet running gunicorn/flask
PRIVATE_HOST=<ip address>
PRIVATE_PORT=<port>

# Public server details
# These are the server details included in the email notifying the client of their finished task
# If running locally or with docker will just be localhost 5000
# If running on digital ocean, use the public ip address of the droplet running gunicorn/flask
# If you have set up a domain name, just use that for the public_host
PUBLIC_HOST=<ip address or domain name>
PUBLIC_PORT=<port>

# This is for gunicorn
# Not actually sure why, just make it a unique random string
SECRET_KEY=<random string>
```    

# Docker

To build the docker image, run the following command in the mrnaid/backend directory of the project:
```
docker compose up redis
docker compose up flower
docker compose up worker
docker compose up flask
```

When Starting docker it will read values in your local .env file by default
Unless you either supply your own like this

```
docker compose --env-file .env.custom up worker
```

or suppling variables manually this

```
PRIVATE_HOST=gunicorn docker compose up worker
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