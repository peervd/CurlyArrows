# Quick Start Guide

## For Immediate Testing (Without EntraID)

If you want to test the application quickly without setting up Azure AD:

### Option 1: Mock Authentication Mode

Create a file `app_dev.py`:

```python
from flask import Flask, session, redirect, url_for
from models import db, User
from blueprints.main import main_bp
from blueprints.api import api_bp
import os

app = Flask(__name__, 
            static_folder='static',
            template_folder='templates')

# Basic configuration
app.config['SECRET_KEY'] = 'dev-secret-key'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///curlyarrows.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['OPENAI_API_KEY'] = os.getenv('OPENAI_API_KEY', '')

# Initialize database
db.init_app(app)

# Register blueprints
app.register_blueprint(main_bp)
app.register_blueprint(api_bp, url_prefix='/api')

# Mock authentication
@app.before_request
def mock_auth():
    """Auto-login for development"""
    if 'user' not in session:
        with app.app_context():
            db.create_all()
            
            # Get or create test user
            test_user = User.query.filter_by(email='test@example.com').first()
            if not test_user:
                test_user = User(
                    azure_oid='test-oid-123',
                    email='test@example.com',
                    name='Test Student'
                )
                db.session.add(test_user)
                db.session.commit()
            
            # Set session
            session['user'] = {
                'name': test_user.name,
                'email': test_user.email,
                'oid': test_user.azure_oid
            }
            session['user_id'] = test_user.id

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
```

Then run:
```bash
python app_dev.py
```

### Option 2: Skip Authentication Temporarily

Edit `blueprints/auth.py` and modify the `login_required` decorator:

```python
def login_required(f):
    """Decorator to require authentication"""
    @wraps(f)
    def decorated_function(*args, **kwargs):
        # FOR TESTING ONLY - Remove in production!
        if 'user' not in session:
            from models import User
            test_user = User.query.filter_by(email='test@example.com').first()
            if test_user:
                session['user'] = {
                    'name': test_user.name,
                    'email': test_user.email,
                    'oid': test_user.azure_oid
                }
                session['user_id'] = test_user.id
        return f(*args, **kwargs)
    return decorated_function
```

## Full Setup with EntraID

### 1. Quick Install

```bash
# Run setup script
chmod +x setup.sh
./setup.sh

# Or manually:
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 2. Minimal Configuration

Create `.env`:
```bash
SECRET_KEY=your-secret-key
OPENAI_API_KEY=your-openai-key

# For testing without EntraID:
FLASK_ENV=development
DATABASE_URL=sqlite:///curlyarrows.db

# Add later when ready:
# AZURE_AD_CLIENT_ID=...
# AZURE_AD_CLIENT_SECRET=...
# AZURE_AD_TENANT_ID=...
```

### 3. Run

```bash
python app.py
```

Visit: http://localhost:5000

## Testing Checklist

- [ ] Application starts without errors
- [ ] Can access the sketcher page
- [ ] ChemDoodle canvas loads
- [ ] Can select an exercise
- [ ] Can draw on canvas
- [ ] Can submit for analysis (requires OpenAI key)
- [ ] Submission appears in database
- [ ] Can view profile page

## Quick Database Check

```python
# Python shell
from app import create_app
from models import User, ExerciseSubmission, db

app = create_app()
with app.app_context():
    # Check users
    users = User.query.all()
    print(f"Users: {len(users)}")
    
    # Check submissions
    submissions = ExerciseSubmission.query.all()
    print(f"Submissions: {len(submissions)}")
```

## Common Quick Start Issues

### Port Already in Use
```bash
# Kill process on port 5000
lsof -ti:5000 | xargs kill -9

# Or use different port
python app.py  # Edit app.py to change port
```

### Missing Dependencies
```bash
# Reinstall all dependencies
pip install --upgrade -r requirements.txt
```

### Database Errors
```bash
# Reset database
rm curlyarrows.db
python app.py  # Will recreate
```

### Static Files Not Loading
```bash
# Check structure
ls -R static/

# Should see:
# static/css/
# static/js/
# static/lib/
# static/exe_img/
```

## Next Steps

Once basic functionality works:

1. Set up Azure AD app registration
2. Add real authentication
3. Configure production database
4. Deploy to server
5. Add SSL certificate
6. Configure backup system

## Getting Help

Check these files for more information:
- `README.md` - Full documentation
- `MIGRATION.md` - Migration details
- `config.py` - Configuration options
- `models.py` - Database schema
