from config import REDIS_PASSWORD, REDIS_ADDRESS, RESULT_EXPIRES
broker_url = f"redis://default:{REDIS_PASSWORD}@{REDIS_ADDRESS}:6379/0"
result_backend = f"redis://default:{REDIS_PASSWORD}@{REDIS_ADDRESS}:6379/1"
result_expires = RESULT_EXPIRES