#!/bin/bash

# Non-interactive Azure Container Apps Deployment Script for CurlyArrows
# Uses Azure Cloud Build (no local Docker required)
set -e

# Configuration
RESOURCE_GROUP="eur-rg-curly"
LOCATION="westeurope"
CONTAINER_REGISTRY="euwcrask"
CONTAINER_REGISTRY_RG="euw-rg-ask"
CONTAINER_APP_NAME="curlyarrows-app"
CONTAINER_ENV_NAME="eur-env-curly"
IMAGE_NAME="curlyarrows-flask"
IMAGE_TAG="latest"

# Load environment variables from .env
if [ -f .env ]; then
    export $(cat .env | grep -v '^#' | xargs)
fi

echo "============================================"
echo "CurlyArrows Azure Deployment (Cloud Build)"
echo "============================================"

# Check Azure login
az account show > /dev/null 2>&1 || { echo "Please login: az login"; exit 1; }

# Create resource group
echo "Creating resource group..."
az group create --name $RESOURCE_GROUP --location $LOCATION --output none || true

# Build image in Azure Cloud (no local Docker needed!)
echo "Building Docker image in Azure Cloud..."
az acr build \
    --registry $CONTAINER_REGISTRY \
    --image ${IMAGE_NAME}:${IMAGE_TAG} \
    --platform linux/amd64 \
    .

# Create environment if it doesn't exist
echo "Creating Container Apps environment..."
if ! az containerapp env show --name $CONTAINER_ENV_NAME --resource-group $RESOURCE_GROUP > /dev/null 2>&1; then
    az containerapp env create \
        --name $CONTAINER_ENV_NAME \
        --resource-group $RESOURCE_GROUP \
        --location $LOCATION \
        --output none
fi

# Get ACR credentials
ACR_USERNAME=$(az acr credential show --name $CONTAINER_REGISTRY --query username -o tsv)
ACR_PASSWORD=$(az acr credential show --name $CONTAINER_REGISTRY --query passwords[0].value -o tsv)

# Check if app exists and get actual FQDN
APP_FQDN=$(az containerapp show \
    --name $CONTAINER_APP_NAME \
    --resource-group $RESOURCE_GROUP \
    --query properties.configuration.ingress.fqdn \
    -o tsv 2>/dev/null || echo "")

if [ -n "$APP_FQDN" ]; then
    # App exists - update it
    echo "Updating existing Container App..."
    REDIRECT_URI="https://${APP_FQDN}/auth/callback"

    az containerapp update \
        --name $CONTAINER_APP_NAME \
        --resource-group $RESOURCE_GROUP \
        --image ${CONTAINER_REGISTRY}.azurecr.io/${IMAGE_NAME}:${IMAGE_TAG} \
        --set-env-vars \
            FLASK_ENV=production \
            SECRET_KEY="${SECRET_KEY}" \
            AZURE_AD_CLIENT_ID="${AZURE_AD_CLIENT_ID}" \
            AZURE_AD_CLIENT_SECRET="${AZURE_AD_CLIENT_SECRET}" \
            AZURE_AD_TENANT_ID="${AZURE_AD_TENANT_ID}" \
            AZURE_AD_REDIRECT_URI="${REDIRECT_URI}" \
            OPENAI_API_KEY="${OPENAI_API_KEY}" \
            DATABASE_URL="sqlite:///instance/curlyarrows.db" \
            HOST=0.0.0.0 \
            PORT=8000 \
            DEBUG=False \
        --output none
else
    # App doesn't exist - create it
    echo "Creating new Container App..."

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
            SECRET_KEY="${SECRET_KEY}" \
            AZURE_AD_CLIENT_ID="${AZURE_AD_CLIENT_ID}" \
            AZURE_AD_CLIENT_SECRET="${AZURE_AD_CLIENT_SECRET}" \
            AZURE_AD_TENANT_ID="${AZURE_AD_TENANT_ID}" \
            AZURE_AD_REDIRECT_URI="placeholder" \
            OPENAI_API_KEY="${OPENAI_API_KEY}" \
            DATABASE_URL="sqlite:///instance/curlyarrows.db" \
            HOST=0.0.0.0 \
            PORT=8000 \
            DEBUG=False \
        --output none

    # Get actual FQDN after creation and update redirect URI
    APP_FQDN=$(az containerapp show \
        --name $CONTAINER_APP_NAME \
        --resource-group $RESOURCE_GROUP \
        --query properties.configuration.ingress.fqdn \
        -o tsv)

    REDIRECT_URI="https://${APP_FQDN}/auth/callback"

    echo "Updating redirect URI with actual FQDN..."
    az containerapp update \
        --name $CONTAINER_APP_NAME \
        --resource-group $RESOURCE_GROUP \
        --set-env-vars AZURE_AD_REDIRECT_URI="${REDIRECT_URI}" \
        --output none
fi

# Get final FQDN
APP_FQDN=$(az containerapp show \
    --name $CONTAINER_APP_NAME \
    --resource-group $RESOURCE_GROUP \
    --query properties.configuration.ingress.fqdn \
    -o tsv)

echo ""
echo "============================================"
echo "Deployment Complete!"
echo "============================================"
echo "App URL: https://$APP_FQDN"
echo "Redirect URI: https://$APP_FQDN/auth/callback"
echo ""
echo "Make sure Azure AD redirect URI is set to:"
echo "  https://$APP_FQDN/auth/callback"
echo "============================================"
