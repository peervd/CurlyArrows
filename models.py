"""
Database Models for CurlyArrows Application
"""
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

db = SQLAlchemy()


def init_db(app):
    """Initialize database"""
    db.init_app(app)
    with app.app_context():
        db.create_all()


class User(db.Model):
    """User model for storing user information from EntraID"""
    __tablename__ = 'users'
    
    id = db.Column(db.Integer, primary_key=True)
    azure_oid = db.Column(db.String(255), unique=True, nullable=False, index=True)
    email = db.Column(db.String(255), unique=True, nullable=False)
    name = db.Column(db.String(255), nullable=False)
    is_teacher = db.Column(db.Boolean, default=False)  # Role: teacher or student
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    updated_at = db.Column(db.DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    submissions = db.relationship('ExerciseSubmission',
                                 foreign_keys='ExerciseSubmission.user_id',
                                 backref='user', lazy='dynamic',
                                 cascade='all, delete-orphan')
    reviewed_submissions = db.relationship('ExerciseSubmission',
                                          foreign_keys='ExerciseSubmission.reviewed_by',
                                          backref='reviewer', lazy='dynamic')
    
    def __repr__(self):
        return f'<User {self.email}>'
    
    def to_dict(self):
        return {
            'id': self.id,
            'email': self.email,
            'name': self.name,
            'is_teacher': self.is_teacher,
            'created_at': self.created_at.isoformat(),
        }


class ExerciseSubmission(db.Model):
    """Model for storing student exercise submissions"""
    __tablename__ = 'exercise_submissions'
    
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'), nullable=False, index=True)
    exercise_number = db.Column(db.Integer, nullable=False, index=True)
    
    # Student's answer data
    student_json_code = db.Column(db.Text, nullable=False)  # ChemDoodle JSON
    student_reasoning = db.Column(db.Text)  # Student's explanation
    
    # Analysis results
    feedback = db.Column(db.Text)  # Generated feedback from analysis
    
    # Timestamps
    submitted_at = db.Column(db.DateTime, default=datetime.utcnow, nullable=False, index=True)
    
    # Teacher review (optional)
    teacher_reviewed = db.Column(db.Boolean, default=False)
    teacher_comments = db.Column(db.Text)
    reviewed_at = db.Column(db.DateTime)
    reviewed_by = db.Column(db.Integer, db.ForeignKey('users.id'))
    
    def __repr__(self):
        return f'<ExerciseSubmission {self.id} - Exercise {self.exercise_number} by User {self.user_id}>'
    
    def to_dict(self):
        return {
            'id': self.id,
            'user_id': self.user_id,
            'exercise_number': self.exercise_number,
            'student_reasoning': self.student_reasoning,
            'feedback': self.feedback,
            'submitted_at': self.submitted_at.isoformat(),
            'teacher_reviewed': self.teacher_reviewed,
            'teacher_comments': self.teacher_comments,
            'reviewed_at': self.reviewed_at.isoformat() if self.reviewed_at else None,
        }


class ExerciseMeta(db.Model):
    """Metadata about exercises (optional)"""
    __tablename__ = 'exercise_meta'
    
    id = db.Column(db.Integer, primary_key=True)
    exercise_number = db.Column(db.Integer, unique=True, nullable=False)
    title = db.Column(db.String(255))
    description = db.Column(db.Text)
    difficulty = db.Column(db.String(50))  # easy, medium, hard
    image_path = db.Column(db.String(255))
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    
    def __repr__(self):
        return f'<ExerciseMeta {self.exercise_number}: {self.title}>'
