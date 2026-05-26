"""
Test script to verify logging configuration
Run with: python test_logging.py
"""
import os
import sys

# Set test environment variables
os.environ['FLASK_ENV'] = 'development'
os.environ['LOG_DIR'] = './logs'
os.environ['LOG_LEVEL'] = 'DEBUG'

from app import create_app

def test_logging():
    """Test the logging functionality"""
    print("Creating Flask application...")
    app = create_app('development')

    print("\nTesting log levels...")
    with app.app_context():
        app.logger.debug("This is a DEBUG message")
        app.logger.info("This is an INFO message")
        app.logger.warning("This is a WARNING message")
        app.logger.error("This is an ERROR message")
        app.logger.critical("This is a CRITICAL message")

    print("\nLogging test complete!")
    print(f"Check console output above and log files in: {app.config['LOG_DIR']}")
    print(f"  - {os.path.join(app.config['LOG_DIR'], 'curlyarrows.log')}")
    print(f"  - {os.path.join(app.config['LOG_DIR'], 'curlyarrows_errors.log')}")

if __name__ == '__main__':
    test_logging()
