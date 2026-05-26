#!/bin/bash
# Setup script for CurlyArrows Flask Application

echo "================================================"
echo "CurlyArrows Flask Application - Setup Script"
echo "================================================"
echo ""

# Check Python version
python_version=$(python3 --version 2>&1 | awk '{print $2}')
echo "✓ Python version: $python_version"

# Create virtual environment
echo ""
echo "Creating virtual environment..."
python3 -m venv venv
echo "✓ Virtual environment created"

# Activate virtual environment
echo ""
echo "Activating virtual environment..."
source venv/bin/activate
echo "✓ Virtual environment activated"

# Upgrade pip
echo ""
echo "Upgrading pip..."
pip install --upgrade pip > /dev/null 2>&1
echo "✓ pip upgraded"

# Install dependencies
echo ""
echo "Installing dependencies..."
echo "This may take a few minutes..."
pip install -r requirements.txt
echo "✓ Dependencies installed"

# Copy environment file
echo ""
if [ ! -f .env ]; then
    echo "Creating .env file from template..."
    cp .env.example .env
    echo "✓ .env file created"
    echo ""
    echo "⚠️  IMPORTANT: Edit .env file with your configuration:"
    echo "   - AZURE_AD_CLIENT_ID"
    echo "   - AZURE_AD_CLIENT_SECRET"
    echo "   - AZURE_AD_TENANT_ID"
    echo "   - OPENAI_API_KEY"
else
    echo "✓ .env file already exists"
fi

# Initialize database
echo ""
echo "Initializing database..."
python3 << EOF
from app import create_app
from models import db

app = create_app()
with app.app_context():
    db.create_all()
    print("✓ Database initialized")
EOF

echo ""
echo "================================================"
echo "Setup Complete!"
echo "================================================"
echo ""
echo "Next steps:"
echo "1. Edit .env file with your Azure AD and OpenAI credentials"
echo "2. Run: source venv/bin/activate"
echo "3. Run: python app.py"
echo "4. Open: http://localhost:5000"
echo ""
echo "For production deployment, see README.md"
echo ""
