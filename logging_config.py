"""
Logging configuration for CurlyArrows application
Supports both console and file logging with configurable directory
"""
import logging
import os
from logging.handlers import RotatingFileHandler
from pathlib import Path


def setup_logging(app):
    """
    Configure logging for the Flask application

    Args:
        app: Flask application instance
    """
    # Get configuration
    log_dir = app.config.get('LOG_DIR', '/tmp/logs')
    log_level = app.config.get('LOG_LEVEL', 'INFO')
    max_bytes = app.config.get('LOG_MAX_BYTES', 10 * 1024 * 1024)
    backup_count = app.config.get('LOG_BACKUP_COUNT', 5)

    # Convert string log level to logging constant
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)

    # Create log directory if it doesn't exist
    try:
        Path(log_dir).mkdir(parents=True, exist_ok=True)
        app.logger.info(f"Log directory created/verified: {log_dir}")
    except Exception as e:
        app.logger.warning(f"Could not create log directory {log_dir}: {e}. Falling back to console only.")
        log_dir = None

    # Remove default handlers
    app.logger.handlers.clear()

    # Set log level
    app.logger.setLevel(numeric_level)

    # Create formatters
    detailed_formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s in %(module)s (%(funcName)s:%(lineno)d): %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    console_formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S'
    )

    # Console Handler - always enabled
    console_handler = logging.StreamHandler()
    console_handler.setLevel(numeric_level)
    console_handler.setFormatter(console_formatter)
    app.logger.addHandler(console_handler)

    # File Handlers - only if log directory is available
    if log_dir:
        # Application log file (all messages)
        app_log_file = os.path.join(log_dir, 'curlyarrows.log')
        app_file_handler = RotatingFileHandler(
            app_log_file,
            maxBytes=max_bytes,
            backupCount=backup_count
        )
        app_file_handler.setLevel(numeric_level)
        app_file_handler.setFormatter(detailed_formatter)
        app.logger.addHandler(app_file_handler)

        # Error log file (errors and critical only)
        error_log_file = os.path.join(log_dir, 'curlyarrows_errors.log')
        error_file_handler = RotatingFileHandler(
            error_log_file,
            maxBytes=max_bytes,
            backupCount=backup_count
        )
        error_file_handler.setLevel(logging.ERROR)
        error_file_handler.setFormatter(detailed_formatter)
        app.logger.addHandler(error_file_handler)

        app.logger.info(f"File logging enabled: {app_log_file}")
        app.logger.info(f"Error logging enabled: {error_log_file}")

    # Log startup information
    app.logger.info("="*60)
    app.logger.info(f"CurlyArrows Application Starting")
    app.logger.info(f"Environment: {app.config.get('ENV', 'unknown')}")
    app.logger.info(f"Log Level: {log_level}")
    app.logger.info(f"Debug Mode: {app.debug}")
    if log_dir:
        app.logger.info(f"Log Directory: {log_dir}")
    app.logger.info("="*60)

    # Also configure werkzeug logger (Flask's underlying server)
    werkzeug_logger = logging.getLogger('werkzeug')
    werkzeug_logger.setLevel(numeric_level)

    # Also configure gunicorn loggers if running under gunicorn
    gunicorn_logger = logging.getLogger('gunicorn.error')
    gunicorn_access_logger = logging.getLogger('gunicorn.access')

    if gunicorn_logger.handlers:
        app.logger.handlers.extend(gunicorn_logger.handlers)
        app.logger.setLevel(gunicorn_logger.level)

    return app.logger


def get_logger(name):
    """
    Get a logger instance for a specific module

    Args:
        name: Logger name (usually __name__)

    Returns:
        Logger instance
    """
    return logging.getLogger(name)
