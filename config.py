"""
Configuration for Flask Application
Includes EntraID/Azure AD authentication settings
"""
import os
from datetime import timedelta

class BaseConfig:
    """Base configuration"""
    SECRET_KEY = os.environ.get('SECRET_KEY', 'dev-secret-key-change-in-production')

    # Session configuration
    SESSION_TYPE = 'filesystem'
    SESSION_PERMANENT = False
    PERMANENT_SESSION_LIFETIME = timedelta(hours=24)
    SESSION_COOKIE_SECURE = True  # Set to True with HTTPS
    SESSION_COOKIE_HTTPONLY = True
    SESSION_COOKIE_SAMESITE = 'Lax'  # Lax allows cookies for same-origin requests

    # Database configuration
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL', 'sqlite:///curlyarrows.db')
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    # File upload configuration
    MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB max file size

    # OpenAI configuration (for LLM feedback)
    OPENAI_API_KEY = os.environ.get('OPENAI_API_KEY', '')

    # Logging configuration
    LOG_DIR = os.environ.get('LOG_DIR', '/tmp/logs')
    LOG_LEVEL = os.environ.get('LOG_LEVEL', 'INFO')
    LOG_MAX_BYTES = 10 * 1024 * 1024  # 10MB per log file
    LOG_BACKUP_COUNT = 5  # Keep 5 backup files


class DevelopmentConfig(BaseConfig):
    """Development configuration"""
    DEBUG = True
    HOST = '0.0.0.0'
    PORT = 5123

    # Override session settings for development with self-signed cert
    SESSION_COOKIE_SECURE = False  # Disable for self-signed certificates

    # EntraID/Azure AD configuration (development placeholders)
    AZURE_AD_CLIENT_ID = os.environ.get('AZURE_AD_CLIENT_ID', 'your-client-id')
    AZURE_AD_CLIENT_SECRET = os.environ.get('AZURE_AD_CLIENT_SECRET', 'your-client-secret')
    AZURE_AD_TENANT_ID = os.environ.get('AZURE_AD_TENANT_ID', 'your-tenant-id')
    AZURE_AD_REDIRECT_URI = os.environ.get('AZURE_AD_REDIRECT_URI', 'https://localhost:5123/auth/callback')

    # Azure AD endpoints
    AZURE_AD_AUTHORITY = f"https://login.microsoftonline.com/{os.environ.get('AZURE_AD_TENANT_ID', 'common')}"
    AZURE_AD_SCOPE = ["User.Read"]

    # Development logging - use local directory
    LOG_DIR = os.environ.get('LOG_DIR', './logs')
    LOG_LEVEL = os.environ.get('LOG_LEVEL', 'DEBUG')


class ProductionConfig(BaseConfig):
    """Production configuration"""
    DEBUG = False
    HOST = '0.0.0.0'
    PORT = int(os.environ.get('PORT', 5000))
    
    # EntraID/Azure AD configuration (must be set via environment variables)
    AZURE_AD_CLIENT_ID = os.environ.get('AZURE_AD_CLIENT_ID')
    AZURE_AD_CLIENT_SECRET = os.environ.get('AZURE_AD_CLIENT_SECRET')
    AZURE_AD_TENANT_ID = os.environ.get('AZURE_AD_TENANT_ID')
    AZURE_AD_REDIRECT_URI = os.environ.get('AZURE_AD_REDIRECT_URI')
    
    # Azure AD endpoints
    AZURE_AD_AUTHORITY = f"https://login.microsoftonline.com/{os.environ.get('AZURE_AD_TENANT_ID')}"
    AZURE_AD_SCOPE = ["User.Read"]
    
    # Enforce HTTPS
    SESSION_COOKIE_SECURE = True
    
    # PostgreSQL for production
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL')


class TestingConfig(BaseConfig):
    """Testing configuration"""
    TESTING = True
    SQLALCHEMY_DATABASE_URI = 'sqlite:///:memory:'
    WTF_CSRF_ENABLED = False
    
    # Mock auth for testing
    AZURE_AD_CLIENT_ID = 'test-client-id'
    AZURE_AD_CLIENT_SECRET = 'test-client-secret'
    AZURE_AD_TENANT_ID = 'test-tenant-id'
    AZURE_AD_REDIRECT_URI = 'http://localhost:5000/auth/callback'


config = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
    'testing': TestingConfig,
    'default': DevelopmentConfig
}
