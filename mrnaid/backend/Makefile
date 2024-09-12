ENV_NAME := mRNAid
SHELL := /bin/bash
CONDA_ACTIVATE = eval "$$(conda shell.bash hook)" && conda activate $(ENV_NAME)
REDIS_VERSION := 6.0.9

env-create: 
	conda env update --name $(ENV_NAME) --file environment.yaml

redis-install:
	cd redis; \
	wget https://download.redis.io/releases/redis-$(REDIS_VERSION).tar.gz; \
	tar xzf redis-$(REDIS_VERSION).tar.gz; \
	cd redis-$(REDIS_VERSION); \
	make;
	cd ../../

redis-run:
	redis/redis-$(REDIS_VERSION)/src/redis-server --bind 0.0.0.0

uwsgi-run:
	$(CONDA_ACTIVATE); \
		APP_NAME=backend \
		CELERY_BROKER_URL=redis://127.0.0.1:6379 \
		CELERY_RESULT_BACKEND=redis://127.0.0.1:6379 \
		LOG_FILE=./logs/logs.log \
		BACKEND_OBJECTIVES_DATA=../common/objectives/data \
		uwsgi --ini app.ini:uwsgi-docker

flask-run:
	$(CONDA_ACTIVATE); \
		flask run --host 0.0.0.0

celery-run:
	$(CONDA_ACTIVATE); \
		LOG_FILE=./logs/logs.log \
		BACKEND_OBJECTIVES_DATA=../common/objectives/data \
		celery -A flask_app.tasks worker --loglevel=info

flower-run:
	$(CONDA_ACTIVATE); \
        celery flower -A flask_app.tasks
