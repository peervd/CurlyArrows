# Use Python 3.11 slim image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies required for chemical libraries
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libopenbabel-dev \
    openbabel \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Create instance directory for SQLite database
RUN mkdir -p instance

# Make entrypoint script executable
RUN chmod +x /app/entrypoint.sh

# Expose port 8000 (Azure Container Apps default)
EXPOSE 8000

# Use entrypoint script to initialize DB and start gunicorn
CMD ["/app/entrypoint.sh"]
