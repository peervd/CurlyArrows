# CurlyArrows Azure Deployment Guide

This guide explains how to deploy the CurlyArrows application to Azure Container Apps.

## Prerequisites

1. **Azure CLI** installed and logged in
   ```bash
   az login
   ```

2. **Docker** installed and running

3. **Azure Resources**:
   - Container Registry: `euwcrask` (already exists)
   - Resource Group: `eur-rg-curly` (will be created)
   - Location: `westeurope`

## Deployment Architecture

The application uses:
- **Azure Container Apps** for hosting the containerized Flask application
- **Azure Container Registry** (euwcrask) for storing Docker images
- **SQLite database** (stored in container - consider PostgreSQL for production)
- **Azure AD (EntraID)** for authentication
- **OpenAI API** for mechanism analysis

## Deployment Steps

### 1. Prepare Environment Variables

Before deploying, gather these values:
- `AZURE_AD_CLIENT_ID` - from Azure Portal > App Registrations
- `AZURE_AD_CLIENT_SECRET` - from Azure Portal > App Registrations
- `AZURE_AD_TENANT_ID` - from Azure Portal > App Registrations
- `OPENAI_API_KEY` - from OpenAI platform

### 2. Run Deployment Script

```bash
./deploy-azure.sh
```

The script will:
1. Create the resource group if needed
2. Build the Docker image
3. Push to Azure Container Registry
4. Create Container Apps environment
5. Deploy the container app
6. Configure environment variables

### 3. Update Azure AD Redirect URI

After deployment, update your Azure AD App Registration:

1. Go to Azure Portal > Azure Active Directory > App Registrations
2. Select your application
3. Go to Authentication > Platform configurations > Web
4. Add redirect URI: `https://curlyarrows-app.westeurope.azurecontainerapps.io/auth/callback`
5. Save

### 4. Initialize Database

After first deployment, initialize the database:

```bash
# Get a shell in the running container
az containerapp exec \
    --name curlyarrows-app \
    --resource-group eur-rg-curly \
    --command '/bin/bash'

# Inside the container, run:
python -c 'from app import app, db; app.app_context().push(); db.create_all()'
exit
```

## Manual Deployment Commands

If you prefer manual deployment:

### Build and Push Image

```bash
# Login to ACR
az acr login --name euwcrask

# Build image
docker build -t euwcrask.azurecr.io/curlyarrows-flask:latest .

# Push image
docker push euwcrask.azurecr.io/curlyarrows-flask:latest
```

### Create Container App

```bash
# Create resource group
az group create --name eur-rg-curly --location westeurope

# Create Container Apps environment
az containerapp env create \
    --name eur-env-curly \
    --resource-group eur-rg-curly \
    --location westeurope

# Get ACR credentials
ACR_USERNAME=$(az acr credential show --name euwcrask --query username -o tsv)
ACR_PASSWORD=$(az acr credential show --name euwcrask --query passwords[0].value -o tsv)

# Create container app
az containerapp create \
    --name curlyarrows-app \
    --resource-group eur-rg-curly \
    --environment eur-env-curly \
    --image euwcrask.azurecr.io/curlyarrows-flask:latest \
    --registry-server euwcrask.azurecr.io \
    --registry-username $ACR_USERNAME \
    --registry-password $ACR_PASSWORD \
    --target-port 8000 \
    --ingress external \
    --min-replicas 0 \
    --max-replicas 3 \
    --cpu 1.0 \
    --memory 2Gi \
    --env-vars \
        FLASK_ENV=production \
        SECRET_KEY="your-secret-key" \
        AZURE_AD_CLIENT_ID="your-client-id" \
        AZURE_AD_CLIENT_SECRET="your-client-secret" \
        AZURE_AD_TENANT_ID="your-tenant-id" \
        AZURE_AD_REDIRECT_URI="https://curlyarrows-app.westeurope.azurecontainerapps.io/auth/callback" \
        OPENAI_API_KEY="your-openai-key" \
        DATABASE_URL="sqlite:///instance/curlyarrows.db" \
        HOST=0.0.0.0 \
        PORT=8000 \
        DEBUG=False
```

## Update Existing Deployment

To update an existing deployment:

```bash
# Build and push new image
docker build -t euwcrask.azurecr.io/curlyarrows-flask:latest .
docker push euwcrask.azurecr.io/curlyarrows-flask:latest

# Update container app
az containerapp update \
    --name curlyarrows-app \
    --resource-group eur-rg-curly \
    --image euwcrask.azurecr.io/curlyarrows-flask:latest
```

## Monitoring and Logs

### View application logs
```bash
az containerapp logs show \
    --name curlyarrows-app \
    --resource-group eur-rg-curly \
    --follow
```

### View metrics
```bash
az containerapp show \
    --name curlyarrows-app \
    --resource-group eur-rg-curly \
    --query properties.configuration
```

## Scaling

The app is configured to scale between 0-3 replicas:
- **Min replicas: 0** - Scales to zero when idle (cost savings)
- **Max replicas: 3** - Scales up under load

To modify scaling:
```bash
az containerapp update \
    --name curlyarrows-app \
    --resource-group eur-rg-curly \
    --min-replicas 1 \
    --max-replicas 5
```

## Environment Variables

Update environment variables:
```bash
az containerapp update \
    --name curlyarrows-app \
    --resource-group eur-rg-curly \
    --set-env-vars \
        VARIABLE_NAME=value
```

## Database Considerations

Currently using SQLite which stores data in the container. For production:

**Limitations**:
- Data lost if container restarts
- Not suitable for multiple replicas
- No backup/restore capabilities

**Recommended**: Migrate to Azure Database for PostgreSQL:

1. Create Azure Database for PostgreSQL
2. Update `DATABASE_URL` environment variable
3. Update `requirements.txt` to include `psycopg2-binary`
4. Rebuild and redeploy

## Logging

The application uses a comprehensive logging system that writes to both console and files.

### Log Configuration

Configure logging via environment variables:

```bash
# Log directory (default: /tmp/logs for Azure, ./logs for local)
LOG_DIR=/tmp/logs

# Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default: INFO)
LOG_LEVEL=INFO
```

### Log Files

When file logging is enabled, the application creates:

1. **`curlyarrows.log`** - All application logs (INFO and above)
2. **`curlyarrows_errors.log`** - Errors and critical issues only

Files are automatically rotated at 10MB with 5 backup files retained.

### Viewing Logs in Azure

```bash
# View recent application logs (console output)
az containerapp logs show --name curlyarrows-app --resource-group eur-rg-curly --tail 100

# Follow logs in real-time
az containerapp logs show --name curlyarrows-app --resource-group eur-rg-curly --follow

# View system logs
az containerapp logs show --name curlyarrows-app --resource-group eur-rg-curly --type system
```

### Log Persistence

**Important**: Logs in `/tmp/logs` are ephemeral and will be lost when:
- Container restarts
- New revision is deployed
- Container scales down to 0

For persistent logs in production, consider:
1. **Azure Application Insights** - Recommended for production monitoring
2. **Azure Storage** - Mount a file share for persistent logs
3. **Log streaming service** - Send logs to external service

### Using Logging in Code

```python
from flask import current_app

# In route handlers
current_app.logger.info("User action")
current_app.logger.warning("Warning message")
current_app.logger.error("Error occurred", exc_info=True)

# In other modules
from logging_config import get_logger
logger = get_logger(__name__)
logger.info("Module-specific log")
```

## Troubleshooting

### Container fails to start
```bash
# Check container logs
az containerapp logs show --name curlyarrows-app --resource-group eur-rg-curly --tail 100

# Check revision status
az containerapp revision list --name curlyarrows-app --resource-group eur-rg-curly
```

### Authentication issues
- Verify Azure AD redirect URI is correct
- Check environment variables are set correctly
- Ensure AZURE_AD_REDIRECT_URI matches the actual app URL

### Database issues
- Ensure database is initialized (see step 4)
- Check file permissions in container
- Consider migrating to PostgreSQL for production

## Cost Optimization

- App scales to 0 replicas when idle (no charges for compute)
- Container Registry storage costs
- Minimal network costs for educational use

## Security Notes

1. **Never commit `.env` file** to version control
2. Use Azure Key Vault for production secrets (optional)
3. Regularly rotate `SECRET_KEY` and API keys
4. Monitor Azure AD sign-in logs for suspicious activity
5. Review OpenAI API usage regularly

## Support

For issues:
1. Check application logs
2. Verify Azure AD configuration
3. Test database connectivity
4. Review environment variables
