#!/bin/bash
# MOSAIC-GUI Google Cloud Platform Deployment Script
# Usage: ./deploy-gcp.sh [project-id] [region]

set -e

# Configuration
PROJECT_ID=${1:-$(gcloud config get-value project)}
REGION=${2:-us-central1}
SERVICE_NAME="mosaic-gui"
IMAGE_NAME="gcr.io/${PROJECT_ID}/${SERVICE_NAME}"

echo "=========================================="
echo "MOSAIC-GUI GCP Deployment"
echo "Project ID: ${PROJECT_ID}"
echo "Region: ${REGION}"
echo "Service Name: ${SERVICE_NAME}"
echo "=========================================="

# Check gcloud is configured
if ! gcloud auth list --filter=status:ACTIVE --format="value(account)" &> /dev/null; then
    echo "Error: gcloud CLI is not authenticated. Run 'gcloud auth login' first."
    exit 1
fi

echo ""
echo "Step 1: Enable Required APIs..."
echo "=========================================="

gcloud services enable \
    run.googleapis.com \
    containerregistry.googleapis.com \
    cloudbuild.googleapis.com \
    --project=${PROJECT_ID}

echo ""
echo "Step 2: Building Container Image..."
echo "=========================================="

cd ../..
gcloud builds submit \
    --tag ${IMAGE_NAME}:latest \
    --project=${PROJECT_ID}

echo ""
echo "Step 3: Deploying to Cloud Run..."
echo "=========================================="

gcloud run deploy ${SERVICE_NAME} \
    --image ${IMAGE_NAME}:latest \
    --platform managed \
    --region ${REGION} \
    --allow-unauthenticated \
    --memory 2Gi \
    --cpu 2 \
    --port 8501 \
    --min-instances 0 \
    --max-instances 10 \
    --concurrency 80 \
    --timeout 300 \
    --set-env-vars "MOSAIC_MODE=demo,STREAMLIT_SERVER_PORT=8501,STREAMLIT_SERVER_ADDRESS=0.0.0.0,STREAMLIT_SERVER_HEADLESS=true" \
    --project=${PROJECT_ID}

echo ""
echo "Step 4: Getting Service URL..."
echo "=========================================="

SERVICE_URL=$(gcloud run services describe ${SERVICE_NAME} \
    --platform managed \
    --region ${REGION} \
    --format 'value(status.url)' \
    --project=${PROJECT_ID})

echo ""
echo "=========================================="
echo "Deployment Complete!"
echo "=========================================="
echo ""
echo "Application URL: ${SERVICE_URL}"
echo ""
echo "To set up a custom domain:"
echo "  gcloud run domain-mappings create \\"
echo "    --service ${SERVICE_NAME} \\"
echo "    --domain your-domain.com \\"
echo "    --region ${REGION} \\"
echo "    --project=${PROJECT_ID}"
echo ""
echo "To view logs:"
echo "  gcloud logging read 'resource.type=cloud_run_revision AND resource.labels.service_name=${SERVICE_NAME}' --project=${PROJECT_ID}"
echo ""
echo "To update the service:"
echo "  gcloud run services update ${SERVICE_NAME} --region ${REGION} --project=${PROJECT_ID} --memory 4Gi"
