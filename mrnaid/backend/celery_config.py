from os import getenv

broker_url = getenv('CELERY_BROKER_URL', "redis://redis:6379/0")
result_backend = getenv('CELERY_RESULT_BACKEND', "redis://redis:6379/1")
result_expires = 60*60*24*30 # Wait for 30 days before deleting the result