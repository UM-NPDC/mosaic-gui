#!/bin/bash
# MOSAIC-GUI Deployment Orchestrator
# Interactive script to deploy to various platforms

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_banner() {
    echo -e "${BLUE}"
    echo "╔═══════════════════════════════════════════════════════════════╗"
    echo "║                                                               ║"
    echo "║   ███╗   ███╗ ██████╗ ███████╗ █████╗ ██╗ ██████╗            ║"
    echo "║   ████╗ ████║██╔═══██╗██╔════╝██╔══██╗██║██╔════╝            ║"
    echo "║   ██╔████╔██║██║   ██║███████╗███████║██║██║                 ║"
    echo "║   ██║╚██╔╝██║██║   ██║╚════██║██╔══██║██║██║                 ║"
    echo "║   ██║ ╚═╝ ██║╚██████╔╝███████║██║  ██║██║╚██████╗            ║"
    echo "║   ╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝╚═╝ ╚═════╝            ║"
    echo "║                                                               ║"
    echo "║           Global Deployment Orchestrator v1.0                 ║"
    echo "║                                                               ║"
    echo "╚═══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_menu() {
    echo ""
    echo -e "${GREEN}Select Deployment Platform:${NC}"
    echo ""
    echo "  FREE OPTIONS:"
    echo "    1) Streamlit Community Cloud  - Quick, easy, free"
    echo "    2) Hugging Face Spaces        - ML community, free tier"
    echo ""
    echo "  CLOUD PROVIDERS:"
    echo "    3) AWS (App Runner)           - Serverless, ~\$25/mo"
    echo "    4) AWS (ECS with GPU)         - Full MOSAIC, ~\$200+/mo"
    echo "    5) Google Cloud Run           - Serverless, ~\$20/mo"
    echo "    6) Google GKE with GPU        - Full MOSAIC, ~\$200+/mo"
    echo "    7) Azure Container Apps       - Serverless, ~\$30/mo"
    echo ""
    echo "  BUDGET OPTIONS:"
    echo "    8) DigitalOcean App Platform  - Simple, ~\$12/mo"
    echo "    9) Docker (Self-hosted)       - Your own server"
    echo ""
    echo "  OTHER:"
    echo "    10) View deployment guide"
    echo "    11) Test locally with Docker"
    echo "    0) Exit"
    echo ""
}

check_docker() {
    if ! command -v docker &> /dev/null; then
        echo -e "${RED}Error: Docker is not installed.${NC}"
        echo "Please install Docker first: https://docs.docker.com/get-docker/"
        exit 1
    fi
}

check_git() {
    if ! command -v git &> /dev/null; then
        echo -e "${RED}Error: Git is not installed.${NC}"
        echo "Please install Git first: https://git-scm.com/downloads"
        exit 1
    fi
}

deploy_streamlit_cloud() {
    echo -e "${YELLOW}Deploying to Streamlit Community Cloud...${NC}"
    echo ""
    
    check_git
    
    echo "Steps to complete:"
    echo "1. Push your code to GitHub"
    echo "2. Go to https://share.streamlit.io"
    echo "3. Sign in with GitHub"
    echo "4. Click 'New app'"
    echo "5. Select your repository"
    echo "6. Set main file path to: app.py"
    echo "7. Click 'Deploy'"
    echo ""
    
    read -p "Would you like to initialize a Git repository? (y/n): " init_git
    if [[ $init_git == "y" ]]; then
        if [ ! -d ".git" ]; then
            git init
            git add .
            git commit -m "Initial commit for Streamlit Cloud deployment"
            echo -e "${GREEN}Git repository initialized!${NC}"
        else
            echo "Git repository already exists."
        fi
        
        echo ""
        echo "Next steps:"
        echo "1. Create a new repository on GitHub"
        echo "2. Run: git remote add origin https://github.com/YOUR_USERNAME/mosaic-gui.git"
        echo "3. Run: git push -u origin main"
        echo "4. Visit https://share.streamlit.io to deploy"
    fi
}

deploy_huggingface() {
    echo -e "${YELLOW}Deploying to Hugging Face Spaces...${NC}"
    echo ""
    
    if ! command -v huggingface-cli &> /dev/null; then
        echo "Installing Hugging Face CLI..."
        pip install huggingface_hub
    fi
    
    echo "Steps:"
    echo "1. Login to Hugging Face"
    echo "2. Create a new Space"
    echo "3. Push code to the Space"
    echo ""
    
    read -p "Would you like to login to Hugging Face? (y/n): " hf_login
    if [[ $hf_login == "y" ]]; then
        huggingface-cli login
        
        echo ""
        read -p "Enter your Hugging Face username: " hf_username
        read -p "Enter your Space name (e.g., mosaic-gui): " space_name
        
        # Copy HuggingFace README
        cp deployment/huggingface/README.md ./README_HF.md
        
        echo ""
        echo "To complete deployment:"
        echo "1. Create a new Space at https://huggingface.co/new-space"
        echo "2. Select 'Streamlit' as the SDK"
        echo "3. Clone the space: git clone https://huggingface.co/spaces/${hf_username}/${space_name}"
        echo "4. Copy your files to the space directory"
        echo "5. Push changes: git add . && git commit -m 'Initial commit' && git push"
    fi
}

deploy_aws_apprunner() {
    echo -e "${YELLOW}Deploying to AWS App Runner...${NC}"
    echo ""
    
    if ! command -v aws &> /dev/null; then
        echo -e "${RED}Error: AWS CLI is not installed.${NC}"
        echo "Install from: https://aws.amazon.com/cli/"
        exit 1
    fi
    
    check_docker
    
    read -p "Enter AWS Region (default: us-east-1): " region
    region=${region:-us-east-1}
    
    cd deployment/aws
    chmod +x deploy-aws.sh
    ./deploy-aws.sh production $region
}

deploy_gcp_cloudrun() {
    echo -e "${YELLOW}Deploying to Google Cloud Run...${NC}"
    echo ""
    
    if ! command -v gcloud &> /dev/null; then
        echo -e "${RED}Error: gcloud CLI is not installed.${NC}"
        echo "Install from: https://cloud.google.com/sdk/docs/install"
        exit 1
    fi
    
    read -p "Enter GCP Project ID: " project_id
    read -p "Enter Region (default: us-central1): " region
    region=${region:-us-central1}
    
    cd deployment/gcp
    chmod +x deploy-gcp.sh
    ./deploy-gcp.sh $project_id $region
}

deploy_docker_local() {
    echo -e "${YELLOW}Building and running locally with Docker...${NC}"
    echo ""
    
    check_docker
    
    echo "Building Docker image..."
    docker build -t mosaic-gui:latest .
    
    echo ""
    echo "Starting container..."
    docker run -d \
        --name mosaic-gui \
        -p 8501:8501 \
        -e MOSAIC_MODE=demo \
        mosaic-gui:latest
    
    echo ""
    echo -e "${GREEN}MOSAIC-GUI is now running!${NC}"
    echo ""
    echo "Access the application at: http://localhost:8501"
    echo ""
    echo "Useful commands:"
    echo "  View logs:    docker logs -f mosaic-gui"
    echo "  Stop:         docker stop mosaic-gui"
    echo "  Remove:       docker rm mosaic-gui"
}

deploy_docker_compose() {
    echo -e "${YELLOW}Deploying with Docker Compose...${NC}"
    echo ""
    
    check_docker
    
    if ! command -v docker-compose &> /dev/null; then
        echo "Using 'docker compose' instead of 'docker-compose'..."
        COMPOSE_CMD="docker compose"
    else
        COMPOSE_CMD="docker-compose"
    fi
    
    echo "Select deployment profile:"
    echo "  1) Demo mode (no GPU)"
    echo "  2) Production with Nginx"
    echo "  3) Full with monitoring"
    read -p "Enter choice (1-3): " compose_choice
    
    case $compose_choice in
        1)
            $COMPOSE_CMD up -d mosaic-demo
            ;;
        2)
            $COMPOSE_CMD --profile production up -d
            ;;
        3)
            $COMPOSE_CMD --profile production --profile monitoring up -d
            ;;
        *)
            echo "Invalid choice"
            exit 1
            ;;
    esac
    
    echo ""
    echo -e "${GREEN}Deployment complete!${NC}"
    echo "Access at: http://localhost:8501"
}

view_guide() {
    if command -v less &> /dev/null; then
        less DEPLOYMENT.md
    elif command -v more &> /dev/null; then
        more DEPLOYMENT.md
    else
        cat DEPLOYMENT.md
    fi
}

# Main script
print_banner

while true; do
    print_menu
    read -p "Enter your choice: " choice
    
    case $choice in
        1) deploy_streamlit_cloud ;;
        2) deploy_huggingface ;;
        3) deploy_aws_apprunner ;;
        4) 
            echo "For GPU deployment, use the CloudFormation template:"
            echo "  cd deployment/aws"
            echo "  aws cloudformation deploy --template-file cloudformation-template.yaml --stack-name mosaic-gpu --capabilities CAPABILITY_IAM --parameter-overrides EnableGPU=true InstanceType=g4dn.xlarge"
            ;;
        5) deploy_gcp_cloudrun ;;
        6)
            echo "For GKE GPU deployment:"
            echo "  kubectl apply -f deployment/gcp/kubernetes-deployment.yaml"
            ;;
        7)
            echo "Azure Container Apps deployment:"
            echo "  az containerapp create --name mosaic-gui --resource-group mosaic-rg --image your-registry.azurecr.io/mosaic-gui:latest --target-port 8501 --ingress external"
            ;;
        8)
            echo "DigitalOcean App Platform:"
            echo "  1. Go to https://cloud.digitalocean.com/apps"
            echo "  2. Create App from GitHub"
            echo "  3. Select your repository"
            echo "  4. Configure as Docker application"
            ;;
        9) deploy_docker_compose ;;
        10) view_guide ;;
        11) deploy_docker_local ;;
        0) 
            echo -e "${GREEN}Thank you for using MOSAIC-GUI!${NC}"
            exit 0 
            ;;
        *) 
            echo -e "${RED}Invalid option. Please try again.${NC}"
            ;;
    esac
    
    echo ""
    read -p "Press Enter to continue..."
done
