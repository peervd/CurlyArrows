# Migration Guide: FastAPI to Flask

## Overview

This document describes the conversion from the original FastAPI application to the new Flask application with EntraID authentication.

## Key Changes

### 1. Framework Migration

**Before (FastAPI)**:
```python
from fastapi import FastAPI
app = FastAPI()

@app.get("/")
def read_root():
    return FileResponse("static/index.html")
```

**After (Flask)**:
```python
from flask import Flask, render_template
app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")
```

### 2. Authentication

**Before**: No authentication
**After**: Microsoft EntraID/Azure AD with MSAL

New authentication flow:
1. User clicks "Login"
2. Redirected to Microsoft login page
3. After authentication, redirected back to app
4. User information stored in session
5. All routes require authentication via `@login_required` decorator

### 3. Application Structure

**Before**: Single file (`server.py`)
**After**: Blueprint-based modular architecture

```
app.py                    # Application factory
├── blueprints/
│   ├── auth.py          # Authentication routes
│   ├── main.py          # Main page routes
│   └── api.py           # API endpoints
```

### 4. Templates

**Before**: Static HTML file served directly
**After**: Jinja2 templates with Bootstrap 5

Benefits:
- Template inheritance (base.html)
- Dynamic content rendering
- Flash messages for user feedback
- Responsive design with Bootstrap

### 5. API Endpoints

**Before**:
```python
@app.post("/api/analyze")
def api_analyze(p: AnalyzePayload):
    # No authentication
    # No database storage
    return result
```

**After**:
```python
@api_bp.route('/analyze', methods=['POST'])
@login_required
def analyze():
    # User authenticated
    # Store in database
    # Associate with user
    return jsonify(result)
```

### 6. Database Integration

**New Feature**: All submissions are now stored in database

Tables:
- `users`: User information from EntraID
- `exercise_submissions`: Student answers and feedback
- `exercise_meta`: Exercise metadata (optional)

Benefits:
- Track student progress
- Teacher can review all submissions
- Historical data for improvement
- Analytics possibilities

### 7. Session Management

**Before**: Stateless (no sessions)
**After**: Server-side sessions

Configuration:
```python
SESSION_TYPE = 'filesystem'
SESSION_PERMANENT = False
PERMANENT_SESSION_LIFETIME = timedelta(hours=24)
```

### 8. Configuration Management

**Before**: Environment variables only
**After**: Structured configuration classes

```python
config.py
├── BaseConfig
├── DevelopmentConfig
├── ProductionConfig
└── TestingConfig
```

### 9. Static Files

**Before**: Mounted at `/static`
**After**: Same structure, but served via Flask's static file handling

No changes needed to static file paths in HTML.

### 10. JavaScript Updates

**Before**: Direct fetch to `/api/analyze`
**After**: Same endpoint, but includes authentication

The JavaScript automatically handles authentication via session cookies.

## Migration Steps for Existing Data

If you have existing data from the FastAPI version:

1. **No database migration needed** - The old version didn't use a database

2. **Static files** - Copy unchanged:
   ```bash
   cp -r old-app/static/* new-app/static/
   ```

3. **PS module** - Copy unchanged:
   ```bash
   cp -r old-app/ps/* new-app/ps/
   ```

## Testing the Migration

### 1. Verify Static Files
- Check that all exercise images load
- Verify ChemDoodle libraries are present
- Test the sketcher functionality

### 2. Test Authentication
- Login with test account
- Verify redirect flow
- Check session persistence
- Test logout functionality

### 3. Test API Endpoints
- Submit a test mechanism
- Verify database storage
- Check feedback generation
- Review submission history

### 4. Test User Roles
- Create teacher account
- Create student account
- Verify permission differences

## Differences in Behavior

### User Experience

**Before**:
- No login required
- No user tracking
- Results shown but not saved
- No history

**After**:
- Login required via EntraID
- User-specific experience
- All submissions saved
- View submission history
- Teachers can review student work

### Security

**Before**:
- Open access
- No user authentication
- No data persistence

**After**:
- Authenticated access only
- User sessions
- Secure cookie handling
- CSRF protection
- Database-backed storage

### Scalability

**Before**:
- Stateless FastAPI
- No database
- Limited to single session

**After**:
- Database-backed
- Multi-user support
- Teacher/student roles
- Historical data analysis

## Configuration Updates Needed

Update your deployment configuration:

1. **Azure AD App Registration** (new requirement)
   - Register application in Azure Portal
   - Configure redirect URIs
   - Set up client secrets

2. **Database** (new requirement)
   - SQLite for development
   - PostgreSQL recommended for production

3. **Environment Variables** (updated)
   ```bash
   # New required variables
   AZURE_AD_CLIENT_ID=...
   AZURE_AD_CLIENT_SECRET=...
   AZURE_AD_TENANT_ID=...
   DATABASE_URL=...
   
   # Existing variables
   OPENAI_API_KEY=...
   SECRET_KEY=...
   ```

## Rollback Plan

If you need to rollback to FastAPI:

1. Keep the old `server.py` file
2. Restore old `requirements.txt`
3. Run: `uvicorn server:app --reload`

However, note that you will lose:
- Authentication
- Database functionality
- User tracking
- Submission history

## Support for Both Versions

You can run both versions simultaneously on different ports:

**FastAPI** (old):
```bash
cd old-app
uvicorn server:app --port 8000
```

**Flask** (new):
```bash
cd flask-app
python app.py  # Runs on port 5000
```

## Gradual Migration Strategy

1. **Week 1**: Deploy Flask version for testing
2. **Week 2**: Parallel run both versions
3. **Week 3**: Migrate users to Flask version
4. **Week 4**: Deprecate FastAPI version

## Common Issues and Solutions

### Issue: "Invalid redirect URI"
**Solution**: Update Azure AD app registration with correct URI

### Issue: "Database locked"
**Solution**: Use PostgreSQL in production instead of SQLite

### Issue: "Session expired"
**Solution**: Check SESSION_PERMANENT_LIFETIME in config

### Issue: "ChemDoodle not loading"
**Solution**: Verify static files copied correctly

## Performance Considerations

- Flask is synchronous by default (unlike FastAPI's async)
- For production, use Gunicorn with multiple workers
- Consider Redis for session storage in production
- Database connection pooling recommended

## Future Enhancements

Possible additions to the Flask version:

1. **Email Notifications**: Alert teachers of new submissions
2. **Analytics Dashboard**: Track student progress
3. **Batch Processing**: Analyze multiple submissions
4. **Export Functions**: Download submission data
5. **API Documentation**: Swagger/OpenAPI integration
6. **WebSocket Support**: Real-time feedback
7. **Mobile App**: Native iOS/Android clients

## Conclusion

The migration from FastAPI to Flask provides:
- ✓ Enterprise authentication (EntraID)
- ✓ Data persistence
- ✓ User management
- ✓ Teacher oversight
- ✓ Progress tracking
- ✓ Scalable architecture

The trade-off is increased complexity, but with significant benefits for educational use.
