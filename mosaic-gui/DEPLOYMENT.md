# MOSAIC-GUI Global Deployment Guide

This guide covers deploying MOSAIC-GUI globally so anyone can access it. Choose the deployment option that best fits your needs and budget.

## Table of Contents

1. [Deployment Options Overview](#deployment-options-overview)
2. [Option 1: Streamlit Community Cloud (FREE)](#option-1-streamlit-community-cloud-free)
3. [Option 2: Hugging Face Spaces (FREE)](#option-2-hugging-face-spaces-free)
4. [Option 3: AWS Deployment](#option-3-aws-deployment)
5. [Option 4: Google Cloud Platform](#option-4-google-cloud-platform)
6. [Option 5: Azure Deployment](#option-5-azure-deployment)
7. [Option 6: DigitalOcean/Linode (Budget-Friendly)](#option-6-digitaloceanlinode-budget-friendly)
8. [Domain & SSL Setup](#domain--ssl-setup)
9. [Scaling Considerations](#scaling-considerations)
10. [Cost Comparison](#cost-comparison)

---

## Deployment Options Overview

| Platform | Cost | GPU Support | Best For |
|----------|------|-------------|----------|
| Streamlit Cloud | FREE | No | Demo/Testing |
| Hugging Face Spaces | FREE | Limited | Public demos |
| AWS | $50-500+/mo | Yes | Production |
| Google Cloud | $50-500+/mo | Yes | Production |
| Azure | $50-500+/mo | Yes | Enterprise |
| DigitalOcean | $20-100/mo | No | Budget deploy |

---

## Option 1: Streamlit Community Cloud (FREE)

**Best for**: Quick deployment, demos, academic use

### Steps:

1. **Create GitHub Repository**
```bash
# Initialize git in your project
cd mosaic-gui
git init
git add .
git commit -m "Initial commit"

# Create repo on GitHub, then push
git remote add origin https://github.com/YOUR_USERNAME/mosaic-gui.git
git branch -M main
git push -u origin main
```

2. **Deploy to Streamlit Cloud**
   - Go to [share.streamlit.io](https://share.streamlit.io)
   - Sign in with GitHub
   - Click "New app"
   - Select your repository: `YOUR_USERNAME/mosaic-gui`
   - Main file path: `app.py`
   - Click "Deploy"

3. **Configure Secrets** (if needed)
   - In Streamlit Cloud dashboard → Your app → Settings → Secrets
   - Add any API keys or configurations

4. **Custom Domain** (optional)
   - Settings → General → Custom subdomain
   - Set: `mosaic-retrosynthesis` → `mosaic-retrosynthesis.streamlit.app`

### Limitations:
- 1GB RAM limit
- No GPU support
- App sleeps after inactivity
- Demo mode only (no full MOSAIC models)

---

## Option 2: Hugging Face Spaces (FREE)

**Best for**: ML community visibility, easy sharing

### Steps:

1. **Create Hugging Face Account**
   - Go to [huggingface.co](https://huggingface.co)
   - Sign up for free account

2. **Create New Space**
```bash
# Install Hugging Face CLI
pip install huggingface_hub

# Login
huggingface-cli login

# Create space
cd mosaic-gui
```

3. **Create `app.py` for Gradio** (alternative to Streamlit)

Create a new file `app_gradio.py`:
```python
import gradio as gr
# ... (see deployment/huggingface/app_gradio.py)
```

4. **Push to Hugging Face**
```bash
# Create requirements.txt for HF
echo "streamlit" > requirements.txt
echo "rdkit" >> requirements.txt
echo "requests" >> requirements.txt
echo "pandas" >> requirements.txt

# Push to space
git remote add hf https://huggingface.co/spaces/YOUR_USERNAME/mosaic-gui
git push hf main
```

### Limitations:
- 16GB RAM on free tier
- Limited GPU hours on free tier
- Community tier limitations

---

## Option 3: AWS Deployment

**Best for**: Production, full MOSAIC with GPU

### Architecture:
```
Internet → Route 53 → CloudFront → ALB → ECS/EC2 → MOSAIC App
                                              ↓
                                           S3 (models)
```

### A. Using AWS App Runner (Simplest)

1. **Push Docker Image to ECR**
```bash
# Configure AWS CLI
aws configure

# Create ECR repository
aws ecr create-repository --repository-name mosaic-gui

# Login to ECR
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin YOUR_ACCOUNT_ID.dkr.ecr.us-east-1.amazonaws.com

# Build and push
docker build -t mosaic-gui .
docker tag mosaic-gui:latest YOUR_ACCOUNT_ID.dkr.ecr.us-east-1.amazonaws.com/mosaic-gui:latest
docker push YOUR_ACCOUNT_ID.dkr.ecr.us-east-1.amazonaws.com/mosaic-gui:latest
```

2. **Create App Runner Service**
```bash
aws apprunner create-service \
  --service-name mosaic-gui \
  --source-configuration '{
    "ImageRepository": {
      "ImageIdentifier": "YOUR_ACCOUNT_ID.dkr.ecr.us-east-1.amazonaws.com/mosaic-gui:latest",
      "ImageRepositoryType": "ECR"
    }
  }' \
  --instance-configuration '{
    "Cpu": "1024",
    "Memory": "2048"
  }'
```

### B. Using ECS with GPU (Full MOSAIC)

See `deployment/aws/ecs-task-definition.json` for complete configuration.

### Cost Estimate:
- App Runner: ~$25-50/month
- ECS with GPU: ~$200-500/month
- S3 storage: ~$5-20/month

---

## Option 4: Google Cloud Platform

**Best for**: Tight integration with Google services, Vertex AI

### A. Using Cloud Run (Serverless)

1. **Enable APIs**
```bash
gcloud services enable run.googleapis.com
gcloud services enable containerregistry.googleapis.com
```

2. **Build and Deploy**
```bash
# Build container
gcloud builds submit --tag gcr.io/YOUR_PROJECT_ID/mosaic-gui

# Deploy to Cloud Run
gcloud run deploy mosaic-gui \
  --image gcr.io/YOUR_PROJECT_ID/mosaic-gui \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --memory 2Gi \
  --cpu 2 \
  --port 8501
```

3. **Custom Domain**
```bash
gcloud run domain-mappings create \
  --service mosaic-gui \
  --domain mosaic.yourdomain.com \
  --region us-central1
```

### B. Using GKE with GPU

See `deployment/gcp/kubernetes-deployment.yaml`

### Cost Estimate:
- Cloud Run: ~$20-50/month
- GKE with GPU: ~$200-400/month

---

## Option 5: Azure Deployment

**Best for**: Enterprise, Microsoft ecosystem

### Using Azure Container Apps

1. **Create Resource Group**
```bash
az group create --name mosaic-rg --location eastus
```

2. **Create Container Registry**
```bash
az acr create --resource-group mosaic-rg --name mosaicregistry --sku Basic
az acr login --name mosaicregistry

# Build and push
az acr build --registry mosaicregistry --image mosaic-gui:latest .
```

3. **Deploy Container App**
```bash
az containerapp create \
  --name mosaic-gui \
  --resource-group mosaic-rg \
  --image mosaicregistry.azurecr.io/mosaic-gui:latest \
  --target-port 8501 \
  --ingress external \
  --cpu 1 \
  --memory 2Gi
```

---

## Option 6: DigitalOcean/Linode (Budget-Friendly)

**Best for**: Cost-effective production deployment

### DigitalOcean App Platform

1. **Connect GitHub Repository**
   - Go to [cloud.digitalocean.com](https://cloud.digitalocean.com)
   - Apps → Create App → GitHub
   - Select your repository

2. **Configure App**
```yaml
# app.yaml (in repo root)
name: mosaic-gui
services:
  - name: web
    github:
      repo: YOUR_USERNAME/mosaic-gui
      branch: main
    dockerfile_path: Dockerfile
    http_port: 8501
    instance_size_slug: professional-xs
    instance_count: 1
```

3. **Deploy**
```bash
doctl apps create --spec app.yaml
```

### Cost Estimate:
- Basic: $12/month
- Professional: $25/month

---

## Domain & SSL Setup

### 1. Purchase Domain
- Namecheap, GoDaddy, Google Domains, Cloudflare

### 2. Configure DNS
```
Type: A
Name: mosaic (or @)
Value: YOUR_SERVER_IP

Type: CNAME
Name: www
Value: mosaic.yourdomain.com
```

### 3. SSL with Let's Encrypt
```bash
# Install certbot
sudo apt install certbot python3-certbot-nginx

# Get certificate
sudo certbot --nginx -d mosaic.yourdomain.com -d www.mosaic.yourdomain.com

# Auto-renewal
sudo certbot renew --dry-run
```

---

## Scaling Considerations

### Horizontal Scaling
```yaml
# docker-compose.scale.yml
services:
  mosaic:
    deploy:
      replicas: 3
      resources:
        limits:
          cpus: '2'
          memory: 4G
```

### Load Balancing
- Use cloud provider's load balancer
- Or deploy Traefik/HAProxy

### Caching
- Redis for session management
- CDN for static assets

### Database (if needed)
- PostgreSQL for user data
- MongoDB for predictions history

---

## Cost Comparison

| Platform | Demo Mode | Full MOSAIC (GPU) |
|----------|-----------|-------------------|
| Streamlit Cloud | FREE | N/A |
| Hugging Face | FREE | ~$50/mo |
| AWS | $25-50/mo | $200-500/mo |
| GCP | $20-50/mo | $200-400/mo |
| Azure | $30-60/mo | $250-500/mo |
| DigitalOcean | $12-25/mo | N/A |

---

## Quick Start Recommendation

### For Academic/Demo Use:
1. Deploy to **Streamlit Community Cloud** (free)
2. Share URL with users

### For Production Use:
1. Deploy to **AWS App Runner** or **Google Cloud Run**
2. Add custom domain with SSL
3. Set up monitoring

### For Full MOSAIC with GPU:
1. Deploy to **AWS ECS** with GPU instances
2. Store models on S3
3. Use spot instances to reduce costs

---

## Support & Monitoring

### Set Up Monitoring
```bash
# Add to docker-compose
docker-compose --profile monitoring up -d
```

### Access Dashboards
- Grafana: http://your-server:3000
- Prometheus: http://your-server:9090

### Alerts
Configure in Grafana for:
- High memory usage
- Error rates
- Response times

---

## Need Help?

- Open an issue on GitHub
- Contact: your-email@domain.com
- Documentation: https://mosaic-docs.yourdomain.com
