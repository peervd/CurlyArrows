"""
API Blueprint
Handles the analysis API endpoint and stores student results
"""
from flask import Blueprint, request, jsonify, session, current_app
from blueprints.auth import login_required
from datetime import datetime
import os

api_bp = Blueprint('api', __name__)


@api_bp.route('/analyze', methods=['POST'])
@login_required
def analyze():
    """
    Analyze student's mechanism drawing
    Stores the result in database for teacher review and LLM feedback
    """
    try:
        # Get request data
        data = request.get_json()
        
        if not data:
            return jsonify({'error': 'No data provided'}), 400
        
        exercise = data.get('exercise')
        student_json_code = data.get('student_json_code')
        student_reasoning = data.get('student_reasoning', '')

        # Debug: Log received data
        current_app.logger.info(f"[DEBUG API] Received exercise={exercise}, type={type(exercise)}")

        # Validate required fields
        if not exercise or not student_json_code:
            return jsonify({'error': 'Missing required fields: exercise and student_json_code'}), 400
        
        # Get user info from session
        user_id = session.get('user_id')

        # Debug: Log session data
        current_app.logger.info(f"API analyze called: session keys={list(session.keys())}, user_id={user_id}")

        if not user_id:
            return jsonify({'error': 'User not authenticated'}), 401
        
        # Import the analysis function from ps module
        from ps.operate_analysis import analyze as analyze_orc
        
        # Get OpenAI key from config
        openai_key = data.get('openai_key') or current_app.config.get('OPENAI_API_KEY') or False
        
        # Perform analysis
        result = analyze_orc(
            openai_key=openai_key,
            exersice=exercise,  # Note: function uses 'exersice' parameter name
            student_json_code=student_json_code,
            student_reasoning=student_reasoning or False,
        )
        
        # Store result in database
        from models import ExerciseSubmission, db
        
        submission = ExerciseSubmission(
            user_id=user_id,
            exercise_number=exercise,
            student_json_code=student_json_code,
            student_reasoning=student_reasoning,
            feedback=result if isinstance(result, str) else str(result),
            submitted_at=datetime.utcnow()
        )
        
        db.session.add(submission)
        db.session.commit()
        
        # Return result
        if isinstance(result, str):
            return jsonify({
                'feedback': result,
                'submission_id': submission.id
            })
        
        return jsonify({
            'result': result,
            'submission_id': submission.id
        })
        
    except Exception as e:
        import traceback
        current_app.logger.error(f"Analysis error: {str(e)}")
        current_app.logger.error(f"Full traceback:\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 400


@api_bp.route('/submissions', methods=['GET'])
@login_required
def get_submissions():
    """
    Get user's submission history
    Teachers can see all submissions, students only their own
    """
    try:
        from models import ExerciseSubmission, User
        
        user_id = session.get('user_id')
        user = User.query.get(user_id)
        
        # Check if user is a teacher (you'll need to implement role checking)
        is_teacher = user.is_teacher if user else False
        
        if is_teacher:
            # Teachers can see all submissions
            student_id = request.args.get('student_id')
            if student_id:
                submissions = ExerciseSubmission.query.filter_by(user_id=student_id).order_by(
                    ExerciseSubmission.submitted_at.desc()
                ).all()
            else:
                submissions = ExerciseSubmission.query.order_by(
                    ExerciseSubmission.submitted_at.desc()
                ).all()
        else:
            # Students only see their own submissions
            submissions = ExerciseSubmission.query.filter_by(user_id=user_id).order_by(
                ExerciseSubmission.submitted_at.desc()
            ).all()
        
        return jsonify({
            'submissions': [
                {
                    'id': sub.id,
                    'exercise_number': sub.exercise_number,
                    'submitted_at': sub.submitted_at.isoformat(),
                    'feedback': sub.feedback,
                    'student_name': sub.user.name if is_teacher else None
                }
                for sub in submissions
            ]
        })
        
    except Exception as e:
        current_app.logger.error(f"Error fetching submissions: {str(e)}")
        return jsonify({'error': str(e)}), 400


@api_bp.route('/submission/<int:submission_id>', methods=['GET'])
@login_required
def get_submission(submission_id):
    """Get details of a specific submission"""
    try:
        from models import ExerciseSubmission, User
        
        user_id = session.get('user_id')
        user = User.query.get(user_id)
        
        submission = ExerciseSubmission.query.get_or_404(submission_id)
        
        # Check permissions
        if not user.is_teacher and submission.user_id != user_id:
            return jsonify({'error': 'Unauthorized'}), 403
        
        return jsonify({
            'id': submission.id,
            'exercise_number': submission.exercise_number,
            'student_json_code': submission.student_json_code,
            'student_reasoning': submission.student_reasoning,
            'feedback': submission.feedback,
            'submitted_at': submission.submitted_at.isoformat(),
            'student_name': submission.user.name
        })
        
    except Exception as e:
        current_app.logger.error(f"Error fetching submission: {str(e)}")
        return jsonify({'error': str(e)}), 400
