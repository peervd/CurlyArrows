#!/bin/bash
set -e

# Use /tmp for SQLite database (writable in Azure Container Apps)
export DATABASE_URL="sqlite:////tmp/curlyarrows.db"

# Initialize database if it doesn't exist
if [ ! -f /tmp/curlyarrows.db ]; then
    echo "Initializing database..."
    python -c 'from app import create_app; from models import db; app = create_app("production"); app.app_context().push(); db.create_all()'
    echo "Database initialized successfully"
else
    echo "Database already exists"
fi

# Start gunicorn with wsgi module
exec gunicorn --bind 0.0.0.0:8000 --workers 2 --timeout 120 --access-logfile - --error-logfile - wsgi:app
