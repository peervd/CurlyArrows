# Implementation Checklist

Use this checklist to implement and deploy the CurlyArrows Flask application.

## Phase 1: Initial Setup (Day 1)

### Local Development Environment

- [ ] Extract the archive: `tar -xzf curlyarrows-flask.tar.gz`
- [ ] Navigate to directory: `cd curlyarrows-flask`
- [ ] Review README.md and QUICKSTART.md
- [ ] Create Python virtual environment: `python3 -m venv venv`
- [ ] Activate virtual environment: `source venv/bin/activate`
- [ ] Install dependencies: `pip install -r requirements.txt`
- [ ] Copy environment template: `cp .env.example .env`
- [ ] Generate SECRET_KEY: `python -c "import secrets; print(secrets.token_hex(32))"`
- [ ] Update SECRET_KEY in .env file

### Initial Test (Without EntraID)

- [ ] Set FLASK_ENV=development in .env
- [ ] Set DATABASE_URL=sqlite:///curlyarrows.db in .env
- [ ] Run: `python app.py`
- [ ] Access http://localhost:5000
- [ ] Verify static files load (ChemDoodle, images)
- [ ] Check console for errors

## Phase 2: Azure AD Setup (Days 2-3)

### Azure Portal Configuration

- [ ] Log in to Azure Portal (https://portal.azure.com)
- [ ] Navigate to Azure Active Directory
- [ ] Go to "App registrations"
- [ ] Click "New registration"

### App Registration Details

- [ ] Set Name: "CurlyArrows Educational Platform"
- [ ] Choose: "Accounts in this organizational directory only"
- [ ] Set Redirect URI: 
  - Type: Web
  - URL: `http://localhost:5000/auth/callback` (development)
- [ ] Click "Register"

### Collect Credentials

- [ ] Copy Application (client) ID → AZURE_AD_CLIENT_ID
- [ ] Copy Directory (tenant) ID → AZURE_AD_TENANT_ID
- [ ] Go to "Certificates & secrets"
- [ ] Click "New client secret"
- [ ] Add description: "CurlyArrows App Secret"
- [ ] Set expiration (recommend: 24 months)
- [ ] Copy the Value (not Secret ID) → AZURE_AD_CLIENT_SECRET
- [ ] ⚠️ Save this immediately - it won't be shown again!

### API Permissions

- [ ] Go to "API permissions"
- [ ] Click "Add a permission"
- [ ] Select "Microsoft Graph"
- [ ] Select "Delegated permissions"
- [ ] Add: "User.Read"
- [ ] Click "Grant admin consent" (if you have permission)

### Update .env File

- [ ] Set AZURE_AD_CLIENT_ID=<your-client-id>
- [ ] Set AZURE_AD_CLIENT_SECRET=<your-secret>
- [ ] Set AZURE_AD_TENANT_ID=<your-tenant-id>
- [ ] Set AZURE_AD_REDIRECT_URI=http://localhost:5000/auth/callback

## Phase 3: OpenAI Integration (Day 3)

### Get OpenAI API Key

- [ ] Go to https://platform.openai.com
- [ ] Sign in or create account
- [ ] Navigate to API keys
- [ ] Create new secret key
- [ ] Copy the key → OPENAI_API_KEY
- [ ] Set in .env file

### Test AI Feedback

- [ ] Restart Flask application
- [ ] Login with Azure AD credentials
- [ ] Select an exercise
- [ ] Draw a simple mechanism
- [ ] Add reasoning text
- [ ] Click "ANALYZE"
- [ ] Verify feedback is generated

## Phase 4: Database Setup (Day 4)

### Development (SQLite)

- [ ] Verify DATABASE_URL=sqlite:///curlyarrows.db in .env
- [ ] Run app to auto-create database
- [ ] Check that curlyarrows.db file is created
- [ ] Verify tables exist

### Check Database

```python
from app import create_app
from models import db, User, ExerciseSubmission

app = create_app()
with app.app_context():
    # List all tables
    print(db.metadata.tables.keys())
    
    # Count users
    user_count = User.query.count()
    print(f"Users: {user_count}")
    
    # Count submissions
    submission_count = ExerciseSubmission.query.count()
    print(f"Submissions: {submission_count}")
```

- [ ] Run database check
- [ ] Verify tables created correctly

## Phase 5: User Management (Day 5)

### Create Test Users

- [ ] Login with your Azure AD account
- [ ] Verify user is created in database
- [ ] Have colleagues login to create more users

### Set Teacher Roles

```python
from app import create_app
from models import db, User

app = create_app()
with app.app_context():
    # Make user a teacher
    user = User.query.filter_by(email='teacher@example.com').first()
    if user:
        user.is_teacher = True
        db.session.commit()
        print(f"Set {user.name} as teacher")
```

- [ ] Identify teacher accounts
- [ ] Set is_teacher flag for each
- [ ] Verify teachers can see all submissions

## Phase 6: Testing (Days 6-7)

### Functional Testing

- [ ] Test login flow (new user)
- [ ] Test login flow (existing user)
- [ ] Test exercise selection
- [ ] Test mechanism drawing
- [ ] Test submission
- [ ] Test feedback display
- [ ] Test submission history
- [ ] Test profile page
- [ ] Test logout

### Role Testing

As Student:
- [ ] Can only see own submissions
- [ ] Cannot access other users' data
- [ ] Can submit exercises
- [ ] Can view own history

As Teacher:
- [ ] Can see all submissions
- [ ] Can filter by student
- [ ] Can add comments
- [ ] Can access analytics

### Error Testing

- [ ] Test with invalid Azure AD credentials
- [ ] Test without OpenAI key
- [ ] Test with malformed mechanism data
- [ ] Test with very long reasoning text
- [ ] Test session timeout
- [ ] Test concurrent users

## Phase 7: Production Preparation (Week 2)

### Database Migration

- [ ] Set up PostgreSQL server
- [ ] Create database: `createdb curlyarrows`
- [ ] Update DATABASE_URL in .env (production)
- [ ] Test connection
- [ ] Migrate data (if any)

### Update Azure AD Redirect URI

- [ ] Go to Azure Portal
- [ ] Navigate to your App registration
- [ ] Go to "Authentication"
- [ ] Add production redirect URI: `https://your-domain.com/auth/callback`
- [ ] Save changes
- [ ] Update AZURE_AD_REDIRECT_URI in production .env

### Security Hardening

- [ ] Change SECRET_KEY for production
- [ ] Set FLASK_ENV=production
- [ ] Enable HTTPS
- [ ] Set SESSION_COOKIE_SECURE=True
- [ ] Configure session timeout
- [ ] Set up Redis for sessions (recommended)
- [ ] Enable HSTS headers
- [ ] Configure CORS properly
- [ ] Set up rate limiting

### Server Setup

- [ ] Install Nginx
- [ ] Configure SSL certificate (Let's Encrypt)
- [ ] Install Gunicorn
- [ ] Create systemd service file
- [ ] Configure Nginx as reverse proxy
- [ ] Set up log rotation
- [ ] Configure firewall

### Deployment

- [ ] Clone repository to server
- [ ] Create virtual environment
- [ ] Install dependencies
- [ ] Configure .env for production
- [ ] Test with Gunicorn: `gunicorn -w 4 -b 0.0.0.0:8000 app:app`
- [ ] Start systemd service
- [ ] Verify HTTPS access
- [ ] Test all functionality

## Phase 8: Monitoring & Maintenance (Ongoing)

### Monitoring Setup

- [ ] Set up application logging
- [ ] Configure error tracking (e.g., Sentry)
- [ ] Monitor disk space (database growth)
- [ ] Monitor memory usage
- [ ] Set up uptime monitoring
- [ ] Configure email alerts

### Backup Strategy

- [ ] Set up daily database backups
- [ ] Test restore procedure
- [ ] Configure off-site backup storage
- [ ] Document recovery process
- [ ] Set retention policy

### Regular Maintenance

Weekly:
- [ ] Review application logs
- [ ] Check error rates
- [ ] Monitor database size
- [ ] Review user feedback

Monthly:
- [ ] Update Python dependencies
- [ ] Review security advisories
- [ ] Check SSL certificate expiry
- [ ] Analyze usage statistics
- [ ] Review and rotate API keys

Quarterly:
- [ ] Test backup restoration
- [ ] Review Azure AD token expiry
- [ ] Update documentation
- [ ] Conduct security audit

## Phase 9: User Training (Week 3)

### Documentation

- [ ] Create user manual for students
- [ ] Create teacher guide
- [ ] Record tutorial videos
- [ ] Prepare FAQ document
- [ ] Create troubleshooting guide

### Training Sessions

- [ ] Schedule student training session
- [ ] Schedule teacher training session
- [ ] Prepare demo exercises
- [ ] Create training materials
- [ ] Collect feedback

## Phase 10: Go Live (Week 4)

### Pre-Launch

- [ ] Final security review
- [ ] Performance testing
- [ ] Load testing (if expecting many users)
- [ ] Backup verification
- [ ] Communication to users

### Launch Day

- [ ] Announce to students and teachers
- [ ] Monitor for issues
- [ ] Be available for support
- [ ] Collect initial feedback
- [ ] Document any issues

### Post-Launch

- [ ] Review first week logs
- [ ] Address any bugs
- [ ] Collect user feedback
- [ ] Plan improvements
- [ ] Update documentation

## Success Criteria

The implementation is successful when:

- [ ] Users can login with institutional credentials
- [ ] Students can complete and submit exercises
- [ ] AI feedback is generated correctly
- [ ] Teachers can review submissions
- [ ] Data is persisted in database
- [ ] Application is stable and responsive
- [ ] Security measures are in place
- [ ] Backup and recovery tested
- [ ] Users are trained
- [ ] Documentation is complete

## Rollback Plan

If issues arise:

1. [ ] Keep old FastAPI version available
2. [ ] Document rollback procedure
3. [ ] Have database snapshot before migration
4. [ ] Test rollback procedure
5. [ ] Communicate clearly with users

## Support Contacts

Document key contacts:

- [ ] Azure AD administrator
- [ ] Database administrator
- [ ] Network/security team
- [ ] Application developers
- [ ] OpenAI support
- [ ] End user support contact

## Notes

Use this space to track progress and issues:

```
Date: _____________
Progress: _________________________________________________________
Issues: ___________________________________________________________
Next steps: _______________________________________________________
```

---

**Remember**: Take backups before each major change!
