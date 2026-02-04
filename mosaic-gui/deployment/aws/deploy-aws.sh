#!/bin/bash
# MOSAIC-GUI AWS Deployment Script
# Usage: ./deploy-aws.sh [environment] [region]

set -e

# Configuration
ENVIRONMENT=${1:-production}
REGION=${2:-us-east-1}
STACK_NAME="mosaic-gui-${ENVIRONMENT}"
ECR_REPO_NAME="${STACK_NAME}-mosaic-gui"

echo "=========================================="
echo "MOSAIC-GUI AWS Deployment"
echo "Environment: ${ENVIRONMENT}"
echo "Region: ${REGION}"
echo "Stack Name: ${STACK_NAME}"
echo "=========================================="

# Check AWS CLI is configured
if ! aws sts get-caller-identity &> /dev/null; then
    echo "Error: AWS CLI is not configured. Run 'aws configure' first."
    exit 1
fi

ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
ECR_URI="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${ECR_REPO_NAME}"

echo ""
echo "Step 1: Deploying CloudFormation Stack..."
echo "=========================================="

aws cloudformation deploy \
    --template-file cloudformation-template.yaml \
    --stack-name ${STACK_NAME} \
    --capabilities CAPABILITY_IAM \
    --parameter-overrides \
        Environment=${ENVIRONMENT} \
    --region ${REGION}

echo ""
echo "Step 2: Getting ECR Repository URI..."
echo "=========================================="

ECR_URI=$(aws cloudformation describe-stacks \
    --stack-name ${STACK_NAME} \
    --query 'Stacks[0].Outputs[?OutputKey==`ECRRepositoryURI`].OutputValue' \
    --output text \
    --region ${REGION})

echo "ECR URI: ${ECR_URI}"

echo ""
echo "Step 3: Logging into ECR..."
echo "=========================================="

aws ecr get-login-password --region ${REGION} | \
    docker login --username AWS --password-stdin ${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com

echo ""
echo "Step 4: Building Docker Image..."
echo "=========================================="

cd ../..
docker build -t mosaic-gui:latest .

echo ""
echo "Step 5: Tagging and Pushing Image..."
echo "=========================================="

docker tag mosaic-gui:latest ${ECR_URI}:latest
docker push ${ECR_URI}:latest

echo ""
echo "Step 6: Updating ECS Service..."
echo "=========================================="

CLUSTER_NAME=$(aws cloudformation describe-stacks \
    --stack-name ${STACK_NAME} \
    --query 'Stacks[0].Outputs[?OutputKey==`ClusterName`].OutputValue' \
    --output text \
    --region ${REGION})

aws ecs update-service \
    --cluster ${CLUSTER_NAME} \
    --service ${STACK_NAME}-service \
    --force-new-deployment \
    --region ${REGION}

echo ""
echo "Step 7: Getting Application URL..."
echo "=========================================="

ALB_DNS=$(aws cloudformation describe-stacks \
    --stack-name ${STACK_NAME} \
    --query 'Stacks[0].Outputs[?OutputKey==`LoadBalancerDNS`].OutputValue' \
    --output text \
    --region ${REGION})

echo ""
echo "=========================================="
echo "Deployment Complete!"
echo "=========================================="
echo ""
echo "Application URL: http://${ALB_DNS}"
echo ""
echo "To set up a custom domain:"
echo "1. Create a Route 53 hosted zone for your domain"
echo "2. Create an A record aliasing to: ${ALB_DNS}"
echo "3. Request an SSL certificate in ACM"
echo "4. Update the ALB listener to use HTTPS"
echo ""
echo "To view logs:"
echo "  aws logs tail /ecs/${STACK_NAME} --follow --region ${REGION}"
echo ""
echo "To scale the service:"
echo "  aws ecs update-service --cluster ${CLUSTER_NAME} --service ${STACK_NAME}-service --desired-count 3 --region ${REGION}"
