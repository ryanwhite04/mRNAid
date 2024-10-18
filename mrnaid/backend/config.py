from dotenv import load_dotenv
from os import getenv

load_dotenv(getenv("ENV", ".env"))

SENDGRID_API_KEY = getenv('SENDGRID_API_KEY')
SENDGRID_EMAIL_USERNAME = getenv('SENDGRID_EMAIL_USERNAME')
NUMBER_OF_ATTEMPTS = 3
REDIS_ADDRESS = getenv('REDIS_ADDRESS', "localhost")
REDIS_PASSWORD = getenv('REDIS_PASSWORD')
PRIVATE_HOST = getenv('PRIVATE_HOST', "localhost")
PRIVATE_PORT = getenv('PRIVATE_PORT', 5000)
PRIVATE_URL = f"http://{PRIVATE_HOST}:{PRIVATE_PORT}"
PUBLIC_HOST = getenv('PUBLIC_HOST', "localhost")
PUBLIC_PORT = getenv('PUBLIC_PORT', 5000)
PUBLIC_URL = f"http://{PUBLIC_HOST}:{PUBLIC_PORT}"# Setting up logger
CELERY_BROKER_URL = getenv('CELERY_BROKER_URL', f"redis://default:{REDIS_PASSWORD}@{REDIS_ADDRESS}:6379/0")
CELERY_RESULT_BACKEND = getenv('CELERY_RESULT_BACKEND', f"redis://default:{REDIS_PASSWORD}@{REDIS_ADDRESS}:6379/1")
RESULT_EXPIRES = 60*60*24*30 # Wait for 30 days before deleting the result
SECRET_KEY = getenv('SECRET_KEY') # No default set because you might forget to change it
LOG_FILE = getenv('LOG_FILE', "./logs/logs.log")