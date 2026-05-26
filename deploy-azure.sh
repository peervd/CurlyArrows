#!/bin/bash

# Azure Container Apps Deployment Script for CurlyArrows
# This script builds and deploys the application to Azure Container Apps

set -e  # Exit on error

# Configuration
RESOURCE_GROUP="eur-rg-curly"
LOCATION="westeurope"  # Change if needed
CONTAINER_REGISTRY="euwcrask"
CONTAINER_REGISTRY_RG="euw-rg-ask"  # Registry is in different resource group
CONTAINER_APP_NAME="curlyarrows-app"
CONTAINER_ENV_NAME="eur-env-curly"
IMAGE_NAME="curlyarrows-flask"
IMAGE_TAG="latest"

echo "============================================"
echo "CurlyArrows Azure Deployment"
echo "============================================"
echo "Resource Group: $RESOURCE_GROUP"
echo "Location: $LOCATION"
echo "Container Registry: $CONTAINER_REGISTRY"
echo "Container App: $CONTAINER_APP_NAME"
echo "============================================"

# Check if logged in to Azure
echo "Checking Azure login status..."
az account show > /dev/null 2>&1 || { echo "Please login to Azure first: az login"; exit 1; }

# Create resource group if it doesn't exist
echo "Creating resource group (if not exists)..."
az group create --name $RESOURCE_GROUP --location $LOCATION --output none || true

# Check if Container Registry exists
echo "Checking container registry..."
if ! az acr show --name $CONTAINER_REGISTRY --resource-group $CONTAINER_REGISTRY_RG > /dev/null 2>&1; then
    echo "ERROR: Container Registry $CONTAINER_REGISTRY does not exist in resource group $CONTAINER_REGISTRY_RG"
    echo "Please verify the registry name and resource group."
    exit 1
fi

# Login to ACR
echo "Logging in to Azure Container Registry..."
az acr login --name $CONTAINER_REGISTRY

# Build and push Docker image for linux/amd64 platform
echo "Building Docker image for linux/amd64..."
docker buildx build --platform linux/amd64 -t ${CONTAINER_REGISTRY}.azurecr.io/${IMAGE_NAME}:${IMAGE_TAG} --push .

# Create Container Apps environment if it doesn't exist
echo "Creating/updating Container Apps environment..."
if ! az containerapp env show --name $CONTAINER_ENV_NAME --resource-group $RESOURCE_GROUP > /dev/null 2>&1; then
    echo "Creating new Container Apps environment..."
    az containerapp env create \
        --name $CONTAINER_ENV_NAME \
        --resource-group $RESOURCE_GROUP \
        --location $LOCATION \
        --output none
else
    echo "Container Apps environment already exists."
fi

# Get ACR credentials
echo "Retrieving ACR credentials..."
ACR_USERNAME=$(az acr credential show --name $CONTAINER_REGISTRY --query username -o tsv)
ACR_PASSWORD=$(az acr credential show --name $CONTAINER_REGISTRY --query passwords[0].value -o tsv)

# Check if Container App exists
if az containerapp show --name $CONTAINER_APP_NAME --resource-group $RESOURCE_GROUP > /dev/null 2>&1; then
    echo "Updating existing Container App..."
    az containerapp update \
        --name $CONTAINER_APP_NAME \
        --resource-group $RESOURCE_GROUP \
        --image ${CONTAINER_REGISTRY}.azurecr.io/${IMAGE_NAME}:${IMAGE_TAG} \
        --output none
else
    echo "Creating new Container App..."

    # Prompt for environment variables
    echo ""
    echo "Please provide the following configuration values:"
    read -p "SECRET_KEY (press Enter for random): " SECRET_KEY
    read -p "AZURE_AD_CLIENT_ID: " AZURE_AD_CLIENT_ID
    read -p "AZURE_AD_CLIENT_SECRET: " AZURE_AD_CLIENT_SECRET
    read -p "AZURE_AD_TENANT_ID: " AZURE_AD_TENANT_ID
    read -p "OPENAI_API_KEY: " OPENAI_API_KEY

    # Generate random secret key if not provided
    if [ -z "$SECRET_KEY" ]; then
        SECRET_KEY=$(python3 -c 'import secrets; print(secrets.token_hex(32))')
        echo "Generated random SECRET_KEY"
    fi

    # Get the app URL (will be available after creation)
    APP_URL="https://${CONTAINER_APP_NAME}.${LOCATION}.azurecontainerapps.io"
    REDIRECT_URI="${APP_URL}/auth/callback"

    echo ""
    echo "IMPORTANT: Update your Azure AD App Registration redirect URI to:"
    echo "  $REDIRECT_URI"
    echo ""
    read -p "Press Enter to continue..."

    az containerapp create \
        --name $CONTAINER_APP_NAME \
        --resource-group $RESOURCE_GROUP \
        --environment $CONTAINER_ENV_NAME \
        --image ${CONTAINER_REGISTRY}.azurecr.io/${IMAGE_NAME}:${IMAGE_TAG} \
        --registry-server ${CONTAINER_REGISTRY}.azurecr.io \
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
            SECRET_KEY="$SECRET_KEY" \
            AZURE_AD_CLIENT_ID="$AZURE_AD_CLIENT_ID" \
            AZURE_AD_CLIENT_SECRET="$AZURE_AD_CLIENT_SECRET" \
            AZURE_AD_TENANT_ID="$AZURE_AD_TENANT_ID" \
            AZURE_AD_REDIRECT_URI="$REDIRECT_URI" \
            OPENAI_API_KEY="$OPENAI_API_KEY" \
            DATABASE_URL="sqlite:///instance/curlyarrows.db" \
            HOST=0.0.0.0 \
            PORT=8000 \
            DEBUG=False \
        --output none
fi

# Get the app URL
APP_FQDN=$(az containerapp show \
    --name $CONTAINER_APP_NAME \
    --resource-group $RESOURCE_GROUP \
    --query properties.configuration.ingress.fqdn \
    -o tsv)

echo ""
echo "============================================"
echo "Deployment Complete!"
echo "============================================"
echo "Application URL: https://$APP_FQDN"
echo ""
echo "IMPORTANT NEXT STEPS:"
echo "1. Update Azure AD App Registration redirect URI to:"
echo "   https://$APP_FQDN/auth/callback"
echo ""
echo "2. Initialize the database by running:"
echo "   az containerapp exec --name $CONTAINER_APP_NAME --resource-group $RESOURCE_GROUP --command '/bin/bash'"
echo "   Then run: python -c 'from app import app, db; app.app_context().push(); db.create_all()'"
echo "============================================"
