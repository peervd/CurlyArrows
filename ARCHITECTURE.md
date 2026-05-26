# CurlyArrows Flask - Architecture Diagram

## System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         Browser (Client)                         │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │          Bootstrap 5 + ChemDoodle Web Components           │ │
│  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐    │ │
│  │  │  Sketcher    │  │  Exercises   │  │   Profile    │    │ │
│  │  │  Interface   │  │     List     │  │     Page     │    │ │
│  │  └──────────────┘  └──────────────┘  └──────────────┘    │ │
│  │                                                             │ │
│  │  JavaScript (base.js, createcanvas-full.js)               │ │
│  └────────────────────────────────────────────────────────────┘ │
└────────────────────────────┬────────────────────────────────────┘
                             │ HTTPS
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                   Flask Application (Server)                     │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │                    app.py (Factory)                        │ │
│  │                                                             │ │
│  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐    │ │
│  │  │   auth_bp    │  │   main_bp    │  │   api_bp     │    │ │
│  │  │              │  │              │  │              │    │ │
│  │  │  /login      │  │  /           │  │  /analyze    │    │ │
│  │  │  /callback   │  │  /exercises  │  │  /submissions│    │ │
│  │  │  /logout     │  │  /instructions│ │              │    │ │
│  │  └──────┬───────┘  └──────────────┘  └──────┬───────┘    │ │
│  │         │                                     │            │ │
│  │         │   ┌──────────────────────────┐    │            │ │
│  │         └───┤  Session Management      ├────┘            │ │
│  │             │  (Flask-Session)         │                 │ │
│  │             └──────────────────────────┘                 │ │
│  └────────────────────────────────────────────────────────────┘ │
└────────────┬──────────────┬─────────────┬──────────────────────┘
             │              │             │
             ▼              ▼             ▼
┌─────────────────┐ ┌──────────────┐ ┌────────────────────┐
│  Microsoft      │ │  Database    │ │  OpenAI API        │
│  EntraID        │ │  (SQLAlchemy)│ │                    │
│  (Azure AD)     │ │              │ │  - GPT Models      │
│                 │ │  Tables:     │ │  - Feedback        │
│  - OAuth2       │ │  - users     │ │  - Analysis        │
│  - User Info    │ │  - exercise_ │ │                    │
│  - Tokens       │ │    submissions│ │                    │
└─────────────────┘ └──────────────┘ └────────────────────┘
```

## Authentication Flow

```
┌──────────┐                                  ┌──────────────┐
│  User    │                                  │  Azure AD    │
└────┬─────┘                                  └──────┬───────┘
     │                                                │
     │  1. Click "Login"                             │
     ├──────────────────────────────────────────────►│
     │                                                │
     │  2. Redirect to Azure AD                      │
     │◄───────────────────────────────────────────────┤
     │                                                │
     │  3. Enter credentials                         │
     ├──────────────────────────────────────────────►│
     │                                                │
     │  4. Authorization code                        │
     │◄───────────────────────────────────────────────┤
     │                                                │
┌────▼─────┐                                  ┌──────▼───────┐
│  Flask   │  5. Exchange code for token     │  Azure AD    │
│  App     ├────────────────────────────────►│              │
│          │                                  │              │
│          │  6. Access token + User info    │              │
│          │◄────────────────────────────────┤              │
└────┬─────┘                                  └──────────────┘
     │
     │  7. Store user in DB
     │  8. Create session
     │  9. Redirect to main page
     │
┌────▼─────┐
│  User    │
│  (Logged │
│   In)    │
└──────────┘
```

## Data Flow for Exercise Submission

```
┌───────────────────────────────────────────────────────────────┐
│  1. Student draws mechanism in ChemDoodle                     │
└───────────────────────────┬───────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────┐
│  2. JavaScript captures JSON representation                   │
└───────────────────────────┬───────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────┐
│  3. POST to /api/analyze with:                                │
│     - exercise_number                                          │
│     - student_json_code                                        │
│     - student_reasoning                                        │
└───────────────────────────┬───────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────┐
│  4. Flask validates authentication (@login_required)          │
└───────────────────────────┬───────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────┐
│  5. Call ps.operate_analysis.analyze()                        │
│     - Parse mechanism                                          │
│     - Compare with model answer                                │
│     - Generate feedback via OpenAI                            │
└───────────────────────────┬───────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────┐
│  6. Store in database:                                        │
│     - Link to user_id                                          │
│     - Save student work                                        │
│     - Save feedback                                            │
│     - Timestamp submission                                     │
└───────────────────────────┬───────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────┐
│  7. Return JSON response to client                            │
└───────────────────────────┬───────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────┐
│  8. Display feedback in UI                                    │
└───────────────────────────────────────────────────────────────┘
```

## Database Schema

```
┌─────────────────────────┐
│       users             │
├─────────────────────────┤
│ id (PK)                 │
│ azure_oid (UNIQUE)      │◄──────┐
│ email (UNIQUE)          │       │
│ name                    │       │
│ is_teacher (BOOLEAN)    │       │
│ created_at              │       │
│ updated_at              │       │
└─────────────────────────┘       │
                                  │
                                  │ Foreign Key
                                  │
┌─────────────────────────────────┴───────────────┐
│       exercise_submissions                      │
├─────────────────────────────────────────────────┤
│ id (PK)                                         │
│ user_id (FK)                                    │
│ exercise_number                                 │
│ student_json_code (TEXT)                        │
│ student_reasoning (TEXT)                        │
│ feedback (TEXT)                                 │
│ submitted_at                                    │
│ teacher_reviewed (BOOLEAN)                      │
│ teacher_comments (TEXT)                         │
│ reviewed_at                                     │
└─────────────────────────────────────────────────┘
```

## Blueprint Organization

```
Flask Application
│
├── auth Blueprint (/auth)
│   ├── /login          → Initiate OAuth flow
│   ├── /callback       → Handle OAuth return
│   ├── /logout         → End session
│   └── /profile        → User profile page
│
├── main Blueprint (/)
│   ├── /               → Main sketcher interface
│   ├── /exercises      → Exercise list
│   └── /instructions   → Help page
│
└── api Blueprint (/api)
    ├── /analyze        → Submit mechanism for analysis
    ├── /submissions    → Get submission history
    └── /submission/:id → Get specific submission
```

## Session Management

```
┌─────────────────────────────────────────────────────────┐
│                    Session Storage                       │
│                   (Flask-Session)                        │
│                                                          │
│  session['user'] = {                                    │
│      'name': 'John Doe',                                │
│      'email': 'john@university.edu',                    │
│      'oid': 'azure-object-id',                          │
│      'given_name': 'John',                              │
│      'surname': 'Doe'                                   │
│  }                                                       │
│                                                          │
│  session['user_id'] = 123  # Database user ID          │
│  session['access_token'] = 'token...'                  │
│                                                          │
│  Lifetime: 24 hours (configurable)                     │
│  Storage: Filesystem (development)                      │
│           Redis (production recommended)                │
└─────────────────────────────────────────────────────────┘
```

## Configuration Hierarchy

```
┌─────────────────────────────────────────┐
│        BaseConfig                        │
│  - SECRET_KEY                           │
│  - SESSION_TYPE                         │
│  - SQLALCHEMY_DATABASE_URI              │
│  - OPENAI_API_KEY                       │
└──────────────┬──────────────────────────┘
               │
      ┌────────┴────────┬─────────────────┐
      │                 │                  │
┌─────▼──────┐   ┌─────▼──────┐   ┌──────▼──────┐
│Development │   │ Production │   │   Testing   │
│            │   │            │   │             │
│- DEBUG: On │   │- DEBUG: Off│   │- In-memory  │
│- SQLite    │   │- PostgreSQL│   │  DB         │
│- HTTP OK   │   │- HTTPS only│   │- Mock auth  │
└────────────┘   └────────────┘   └─────────────┘
```

## Security Layers

```
┌────────────────────────────────────────────────────┐
│  Layer 1: HTTPS/TLS (Production)                  │
└─────────────────┬──────────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────────┐
│  Layer 2: OAuth2 Authentication (EntraID)         │
└─────────────────┬──────────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────────┐
│  Layer 3: Session Management (HTTPOnly, Secure)   │
└─────────────────┬──────────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────────┐
│  Layer 4: CSRF Protection (State tokens)          │
└─────────────────┬──────────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────────┐
│  Layer 5: SQL Injection Prevention (ORM)          │
└─────────────────┬──────────────────────────────────┘
                  │
┌─────────────────▼──────────────────────────────────┐
│  Layer 6: XSS Prevention (Jinja2 auto-escape)     │
└────────────────────────────────────────────────────┘
```

## Deployment Architecture

```
┌─────────────┐
│   Nginx     │  ← Reverse proxy, SSL termination
│   (Port 443)│
└──────┬──────┘
       │
┌──────▼──────┐
│  Gunicorn   │  ← WSGI server (4 workers)
│  (Port 8000)│
└──────┬──────┘
       │
┌──────▼──────┐
│  Flask App  │  ← Application instances
│  (Multiple  │
│   workers)  │
└──────┬──────┘
       │
       ├──────────┬─────────────┬──────────────┐
       │          │             │              │
┌──────▼────┐ ┌──▼──────┐ ┌───▼──────┐ ┌────▼─────┐
│PostgreSQL │ │  Redis  │ │ Azure AD │ │ OpenAI   │
│           │ │(Sessions│ │          │ │   API    │
└───────────┘ └─────────┘ └──────────┘ └──────────┘
```
