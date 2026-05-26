"""
Flask Application Factory for CurlyArrows Educational Application
with Microsoft EntraID Authentication
"""
from flask import Flask
from flask_session import Session
from pathlib import Path
import os


def create_app(config_name='development'):
    """Application factory pattern"""
    app = Flask(__name__,
                static_folder='static',
                template_folder='templates')

    # Load configuration
    app.config.from_object(f'config.{config_name.capitalize()}Config')

    # Setup logging (must be done after config is loaded)
    from logging_config import setup_logging
    setup_logging(app)

    app.logger.info(f"Creating application with config: {config_name}")

    # Initialize Flask-Session for server-side session storage
    Session(app)
    app.logger.info("Flask-Session initialized")

    # Register blueprints
    from blueprints.main import main_bp
    from blueprints.auth import auth_bp
    from blueprints.api import api_bp

    app.register_blueprint(main_bp)
    app.register_blueprint(auth_bp, url_prefix='/auth')
    app.register_blueprint(api_bp, url_prefix='/api')
    app.logger.info("Blueprints registered: main, auth, api")

    # Initialize database
    from models import init_db
    with app.app_context():
        init_db(app)
    app.logger.info("Database initialized")

    return app

if __name__ == '__main__':
    app = create_app(os.getenv('FLASK_ENV', 'development'))

    # SSL context for HTTPS
    ssl_context = ('cert.pem', 'key.pem')

    app.run(
        host=app.config.get('HOST', '0.0.0.0'),
        port=app.config.get('PORT', 5000),
        debug=app.config.get('DEBUG', True),
        ssl_context=ssl_context
    )
