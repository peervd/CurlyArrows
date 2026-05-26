"""
Main Blueprint
Serves the main application pages
"""
from flask import Blueprint, render_template, session, current_app
from blueprints.auth import login_required

main_bp = Blueprint('main', __name__)


@main_bp.route('/')
@login_required
def index():
    """Main page with ChemDoodle sketcher"""
    user = session.get('user', {})
    current_app.logger.info(f"User {user.get('name', 'Unknown')} accessed index page")
    return render_template('index.html', user=user)


@main_bp.route('/instructions')
@login_required
def instructions():
    """Instructions page"""
    return render_template('instructions.html')


@main_bp.route('/exercises')
@login_required
def exercises():
    """List of available exercises"""
    # Load exercise list (could be from database or config)
    exercises = list(range(1, 25))  
    return render_template('exercises.html', exercises=exercises)
