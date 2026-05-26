# CurlyArrows - Flask Educational Platform

Educational application for teaching organic chemistry reaction mechanisms with AI-powered feedback, using Flask and Microsoft EntraID authentication.

## Overview

CurlyArrows is a web-based educational platform that allows students to:
- Draw organic reaction mechanisms using ChemDoodle Web Components
- Receive AI-powered feedback on their mechanisms
- Track their learning progress
- Get personalized instruction from teachers

## Features

- **Microsoft EntraID Authentication**: Secure single sign-on (SSO) with Azure Active Directory
- **ChemDoodle Integration**: Professional chemistry drawing tools
- **AI Feedback**: OpenAI-powered analysis of student mechanisms
- **Database Storage**: All submissions saved for teacher review and progress tracking
- **Blueprint Architecture**: Modular Flask application structure
- **Bootstrap 5 UI**: Modern, responsive interface
- **Teacher Dashboard**: Review student submissions and provide feedback

## Architecture

### Technology Stack

- **Backend**: Flask 3.0
- **Authentication**: Microsoft MSAL (Microsoft Authentication Library)
- **Database**: SQLAlchemy (SQLite for development, PostgreSQL for production)
- **Frontend**: Bootstrap 5, ChemDoodle Web Components
- **AI**: OpenAI API for mechanism analysis
- **Chemistry**: RDKit for molecular analysis

### Application Structure

```
curlyarrows-flask/
├── app.py                  # Application factory
├── config.py               # Configuration settings
├── models.py               # Database models
├── requirements.txt        # Python dependencies
├── .env.example           # Environment variables template
│
├── blueprints/            # Flask blueprints
│   ├── __init__.py
│   ├── auth.py           # EntraID authentication
│   ├── main.py           # Main page routes
│   └── api.py            # Analysis API endpoints
│
├── templates/             # Jinja2 templates
│   ├── base.html         # Base template with Bootstrap
│   ├── index.html        # Main sketcher interface
│   ├── profile.html      # User profile page
│   ├── instructions.html # Instructions page
│   └── exercises.html    # Exercise list
│
├── static/                # Static files
│   ├── css/
│   │   └── styles.css
│   ├── js/
│   │   ├── base.js
│   │   └── createcanvas-full.js
│   ├── lib/              # ChemDoodle libraries
│   └── exe_img/          # Exercise images
│
└── ps/                    # Python modules for analysis
    ├── operate_analysis.py
    ├── mechanism_analysis.py
    └── ...                # Other chemistry analysis modules
```

## Setup Instructions

### 1. Prerequisites

- Python 3.9 or higher
- pip and virtualenv
- Azure AD tenant (for EntraID authentication)
- OpenAI API key (for AI feedback)

### 2. Clone and Install

```bash
# Navigate to the application directory
cd curlyarrows-flask

# Create virtual environment
python3 -m venv venv

# Activate virtual environment
# On macOS/Linux:
source venv/bin/activate
# On Windows:
venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### 3. Configure Environment Variables

Copy the example environment file:

```bash
cp .env.example .env
```

Edit `.env` with your configuration:

```bash
# Flask Configuration
FLASK_ENV=development
SECRET_KEY=your-random-secret-key-here

# Database
DATABASE_URL=sqlite:///curlyarrows.db

# Microsoft EntraID / Azure AD
AZURE_AD_CLIENT_ID=your-client-id
AZURE_AD_CLIENT_SECRET=your-client-secret
AZURE_AD_TENANT_ID=your-tenant-id
AZURE_AD_REDIRECT_URI=http://localhost:5000/auth/callback

# OpenAI
OPENAI_API_KEY=your-openai-api-key
```

### 4. Azure AD App Registration

To enable EntraID authentication, register an application in Azure Portal:

1. Go to [Azure Portal](https://portal.azure.com)
2. Navigate to **Azure Active Directory** > **App registrations**
3. Click **New registration**
4. Configure:
   - **Name**: CurlyArrows
   - **Supported account types**: Accounts in this organizational directory only
   - **Redirect URI**: Web - `http://localhost:5000/auth/callback`

5. After registration:
   - Note the **Application (client) ID** → `AZURE_AD_CLIENT_ID`
   - Note the **Directory (tenant) ID** → `AZURE_AD_TENANT_ID`

6. Create a client secret:
   - Go to **Certificates & secrets**
   - Click **New client secret**
   - Note the **Value** → `AZURE_AD_CLIENT_SECRET`

7. Set API permissions:
   - Go to **API permissions**
   - Add **Microsoft Graph** > **Delegated permissions** > **User.Read**
   - Grant admin consent (if required)

### 5. Initialize Database

```bash
# The database will be created automatically on first run
python app.py
```

Or initialize manually:

```python
from app import create_app
from models import db

app = create_app()
with app.app_context():
    db.create_all()
    print("Database initialized!")
```

### 6. Run the Application

```bash
# Development mode
python app.py

# Or with Flask CLI
flask run
```

The application will be available at: http://localhost:5000

## Configuration

### Development vs Production

The application uses different configurations based on the `FLASK_ENV` variable:

**Development** (`FLASK_ENV=development`):
- Debug mode enabled
- SQLite database
- Detailed error messages
- HTTP allowed for redirect URIs

**Production** (`FLASK_ENV=production`):
- Debug mode disabled
- PostgreSQL database recommended
- HTTPS required
- Enhanced security settings

### Database Options

**SQLite** (Development):
```
DATABASE_URL=sqlite:///curlyarrows.db
```

**PostgreSQL** (Production):
```
DATABASE_URL=postgresql://username:password@localhost/curlyarrows
```

**MySQL**:
```
DATABASE_URL=mysql://username:password@localhost/curlyarrows
```

## Usage

### For Students

1. **Login**: Click "Login" and authenticate with your institutional credentials
2. **Select Exercise**: Choose an exercise from the dropdown menu
3. **Draw Mechanism**: Use the ChemDoodle sketcher to draw your reaction mechanism
4. **Add Reasoning**: Explain your thinking in the reasoning text area
5. **Analyze**: Click "ANALYZE" to receive AI-powered feedback
6. **Review**: Check your profile to see past submissions

### For Teachers

Teachers can:
- View all student submissions via the API endpoint `/api/submissions`
- Access detailed submission data
- Add teacher comments to student work
- Track student progress over time

To designate a user as a teacher, update the database:

```python
from models import User, db
user = User.query.filter_by(email='teacher@example.com').first()
user.is_teacher = True
db.session.commit()
```

## API Endpoints

### Authentication

- `GET /auth/login` - Initiate EntraID login
- `GET /auth/callback` - Handle OAuth2 callback
- `GET /auth/logout` - Logout user
- `GET /auth/profile` - View user profile

### Main Routes

- `GET /` - Main sketcher interface (requires auth)
- `GET /exercises` - List of exercises
- `GET /instructions` - Instructions page

### API Routes

- `POST /api/analyze` - Analyze student mechanism
  ```json
  {
    "exercise": 1,
    "student_json_code": "{ ... }",
    "student_reasoning": "I think..."
  }
  ```

- `GET /api/submissions` - Get user's submissions
- `GET /api/submission/<id>` - Get specific submission details

## Database Schema

### User Table
- `id`: Primary key
- `azure_oid`: Azure AD Object ID (unique)
- `email`: User email
- `name`: Display name
- `is_teacher`: Boolean flag
- `created_at`: Registration timestamp

### ExerciseSubmission Table
- `id`: Primary key
- `user_id`: Foreign key to User
- `exercise_number`: Exercise number (1-14)
- `student_json_code`: ChemDoodle JSON data
- `student_reasoning`: Student's explanation
- `feedback`: AI-generated feedback
- `submitted_at`: Submission timestamp
- `teacher_reviewed`: Boolean flag
- `teacher_comments`: Teacher's feedback
- `reviewed_at`: Review timestamp

## Security Considerations

1. **Session Management**: Server-side sessions using Flask-Session
2. **CSRF Protection**: State tokens in OAuth2 flow
3. **SQL Injection**: Prevented by SQLAlchemy ORM
4. **XSS Protection**: Jinja2 auto-escaping
5. **Secure Cookies**: HTTPOnly and Secure flags in production
6. **HTTPS**: Required in production environment

## Deployment

### Using Gunicorn (Production)

```bash
# Install gunicorn
pip install gunicorn

# Run with gunicorn
gunicorn -w 4 -b 0.0.0.0:5000 app:app
```

### Using Docker

Create a `Dockerfile`:

```dockerfile
FROM python:3.11-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["gunicorn", "-w", "4", "-b", "0.0.0.0:5000", "app:app"]
```

Build and run:

```bash
docker build -t curlyarrows .
docker run -p 5000:5000 --env-file .env curlyarrows
```

### Environment-Specific Settings

Update your production `.env`:

```bash
FLASK_ENV=production
DATABASE_URL=postgresql://user:pass@db-host/curlyarrows
AZURE_AD_REDIRECT_URI=https://your-domain.com/auth/callback
SESSION_COOKIE_SECURE=True
```

## Troubleshooting

### Authentication Issues

**Problem**: "Invalid state parameter" error
- **Solution**: Check that cookies are enabled and session storage is working

**Problem**: "Token acquisition error"
- **Solution**: Verify client ID, secret, and tenant ID are correct

### Database Issues

**Problem**: "No module named 'psycopg2'"
- **Solution**: Install PostgreSQL adapter: `pip install psycopg2-binary`

**Problem**: Database file locked
- **Solution**: Close other connections to the SQLite database

### ChemDoodle Issues

**Problem**: Sketcher not loading
- **Solution**: Verify ChemDoodle library files are in `static/lib/`

## Development

### Adding New Exercises

1. Add exercise image to `static/exe_img/` (e.g., `exercise_15.png`)
2. Update exercise list in templates
3. Add exercise metadata to database (optional)

### Customizing AI Feedback

Edit `ps/operate_analysis.py` to modify the analysis logic and feedback generation.

### Adding New Routes

Create new blueprint files in `blueprints/` directory and register them in `app.py`.

## License

This application is for educational use at Wageningen University & Research.

## Support

For issues or questions, contact the development team or submit an issue in the repository.

## Credits

- **ChemDoodle Web Components**: iChemLabs
- **Flask Framework**: Pallets Projects
- **Bootstrap**: Twitter
- **Microsoft Authentication**: Microsoft
