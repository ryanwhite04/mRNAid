from loguru import logger
import os


def MyLogger(name):
    """ Logger which can be imported to any module """

    logger.remove()

    logger.level("DEBUG", color="<yellow>")
    logger.level("INFO", color="<green>")
    logger.level("ERROR", color="<red>")

    file = os.environ.get('LOG_FILE', "./logs/logs.log")

    if file:
        log_format = "{time:YYYY-MM-DD HH:mm:ss}| {module}| {level}| {message}"
        logger.add(file, format=log_format, level="DEBUG", backtrace=True, rotation="100 MB", compression="zip",
                   enqueue=True)
    logger.enable(name)

    return logger
