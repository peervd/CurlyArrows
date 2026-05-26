"""
Authentication Blueprint
Handles Microsoft EntraID/Azure AD authentication
"""
from flask import Blueprint, redirect, url_for, session, request, current_app, flash
from functools import wraps
import msal
import requests

auth_bp = Blueprint('auth', __name__)


def login_required(f):
    """Decorator to require authentication"""
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if not session.get('user'):
            return redirect(url_for('auth.login', next=request.url))
        return f(*args, **kwargs)
    return decorated_function


def get_msal_app():
    """Create MSAL confidential client application"""
    return msal.ConfidentialClientApplication(
        current_app.config['AZURE_AD_CLIENT_ID'],
        authority=current_app.config['AZURE_AD_AUTHORITY'],
        client_credential=current_app.config['AZURE_AD_CLIENT_SECRET'],
    )


def get_auth_url(state=None):
    """Generate authorization URL for Azure AD"""
    msal_app = get_msal_app()
    auth_url = msal_app.get_authorization_request_url(
        scopes=current_app.config['AZURE_AD_SCOPE'],
        state=state or session.get('state'),
        redirect_uri=current_app.config['AZURE_AD_REDIRECT_URI']
    )
    return auth_url


@auth_bp.route('/login')
def login():
    """Initiate login flow"""
    # Generate state token for CSRF protection
    import uuid
    state = str(uuid.uuid4())
    session['state'] = state
    session['next_url'] = request.args.get('next', url_for('main.index'))
    
    # Redirect to Azure AD for authentication
    auth_url = get_auth_url(state)
    return redirect(auth_url)


@auth_bp.route('/callback')
def callback():
    """Handle callback from Azure AD"""
    current_app.logger.info(f"Callback received: state={request.args.get('state')}, session_state={session.get('state')}")

    # Verify state token
    if request.args.get('state') != session.get('state'):
        current_app.logger.error(f"State mismatch: request={request.args.get('state')}, session={session.get('state')}")
        flash('Invalid state parameter. Please try logging in again.', 'danger')
        return redirect(url_for('main.index'))
    
    # Check for error in response
    if 'error' in request.args:
        flash(f"Authentication error: {request.args.get('error_description', 'Unknown error')}", 'danger')
        return redirect(url_for('main.index'))
    
    # Get authorization code
    code = request.args.get('code')
    if not code:
        flash('No authorization code received', 'danger')
        return redirect(url_for('main.index'))
    
    try:
        # Exchange code for token
        msal_app = get_msal_app()
        result = msal_app.acquire_token_by_authorization_code(
            code,
            scopes=current_app.config['AZURE_AD_SCOPE'],
            redirect_uri=current_app.config['AZURE_AD_REDIRECT_URI']
        )
        
        if 'error' in result:
            flash(f"Token acquisition error: {result.get('error_description', 'Unknown error')}", 'danger')
            return redirect(url_for('main.index'))
        
        # Store user info in session
        session['user'] = {
            'name': result.get('id_token_claims', {}).get('name', 'Unknown'),
            'email': result.get('id_token_claims', {}).get('preferred_username', ''),
            'oid': result.get('id_token_claims', {}).get('oid', ''),  # Object ID in Azure AD
        }
        session['access_token'] = result.get('access_token')
        
        # Get additional user info from Microsoft Graph API (optional)
        try:
            graph_data = requests.get(
                'https://graph.microsoft.com/v1.0/me',
                headers={'Authorization': f"Bearer {result.get('access_token')}"}
            ).json()
            
            session['user']['given_name'] = graph_data.get('givenName', '')
            session['user']['surname'] = graph_data.get('surname', '')
        except:
            pass  # Continue even if Graph API call fails
        
        # Store user in database
        from models import User, db

        current_app.logger.info(f"Looking up user with azure_oid={session['user']['oid']}")

        user = User.query.filter_by(azure_oid=session['user']['oid']).first()
        if not user:
            current_app.logger.info(f"Creating new user: {session['user']['email']}")
            user = User(
                azure_oid=session['user']['oid'],
                email=session['user']['email'],
                name=session['user']['name']
            )
            db.session.add(user)
        else:
            current_app.logger.info(f"Updating existing user: {session['user']['email']}")
            # Update user info if it changed
            user.email = session['user']['email']
            user.name = session['user']['name']

        db.session.commit()
        current_app.logger.info(f"User committed to database: user.id={user.id}")

        session['user_id'] = user.id
        current_app.logger.info(f"User logged in: user_id={user.id}, session keys={list(session.keys())}")

        flash(f"Welcome, {session['user']['name']}!", 'success')

        # Redirect to original destination or home
        next_url = session.pop('next_url', url_for('main.index'))
        return redirect(next_url)
        
    except Exception as e:
        current_app.logger.error(f"Authentication error: {str(e)}")
        flash(f"Authentication failed: {str(e)}", 'danger')
        return redirect(url_for('main.index'))


@auth_bp.route('/logout')
def logout():
    """Logout user"""
    session.clear()
    
    # Redirect to Azure AD logout
    logout_url = (
        f"{current_app.config['AZURE_AD_AUTHORITY']}/oauth2/v2.0/logout"
        f"?post_logout_redirect_uri={url_for('main.index', _external=True)}"
    )
    
    flash('You have been logged out successfully.', 'info')
    return redirect(logout_url)


@auth_bp.route('/profile')
@login_required
def profile():
    """User profile page"""
    from flask import render_template
    return render_template('profile.html', user=session.get('user'))
