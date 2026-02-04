#!/bin/bash
# MOSAIC-GUI Setup Script
# =======================
# This script sets up the MOSAIC-GUI environment

set -e

echo "================================================"
echo "   MOSAIC-GUI Installation Script"
echo "================================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check Python version
echo -n "Checking Python version... "
python_version=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
required_version="3.11"

if [ "$(printf '%s\n' "$required_version" "$python_version" | sort -V | head -n1)" = "$required_version" ]; then
    echo -e "${GREEN}OK ($python_version)${NC}"
else
    echo -e "${YELLOW}Warning: Python $python_version detected. Recommended: 3.11+${NC}"
fi

# Create virtual environment
echo ""
echo "Creating virtual environment..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
    echo -e "${GREEN}Virtual environment created${NC}"
else
    echo -e "${YELLOW}Virtual environment already exists${NC}"
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo ""
echo "Upgrading pip..."
pip install --upgrade pip wheel setuptools

# Install requirements
echo ""
echo "Installing dependencies..."
pip install -r requirements.txt

# Check CUDA availability
echo ""
echo -n "Checking CUDA availability... "
python3 -c "import torch; print('Available' if torch.cuda.is_available() else 'Not available')" 2>/dev/null || echo "PyTorch not yet installed"

# Create directories
echo ""
echo "Creating directory structure..."
mkdir -p Model_Weights
mkdir -p data
mkdir -p outputs
mkdir -p logs

# Download sample data (optional)
echo ""
read -p "Download sample reaction data? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Sample data download would go here..."
    # wget -O data/sample_reactions.csv https://example.com/sample_data.csv
fi

# Check for MOSAIC model weights
echo ""
echo "Checking for MOSAIC model weights..."
if [ -d "Model_Weights/faiss_index" ]; then
    echo -e "${GREEN}FAISS index found${NC}"
else
    echo -e "${YELLOW}FAISS index not found${NC}"
    echo "For full functionality, download model weights from:"
    echo "https://github.com/haoteli/MOSAIC"
fi

# Create config file if not exists
if [ ! -f "config.yaml" ]; then
    echo ""
    echo "Creating default configuration..."
    cat > config.yaml << 'EOF'
# MOSAIC-GUI Configuration
# ========================

# Model paths
model:
  base_path: "./Model_Weights"
  faiss_index: "./Model_Weights/faiss_index"
  kernel_network: "./Model_Weights/kernel_metric_network.pt"
  num_experts: 2498

# Prediction settings
prediction:
  max_steps: 10
  beam_width: 5
  temperature: 0.7
  top_p: 0.9

# Server settings
server:
  host: "localhost"
  port: 8501
  debug: false

# API settings
api:
  pubchem_timeout: 10
  cache_enabled: true
  cache_dir: "./cache"

# Logging
logging:
  level: "INFO"
  file: "./logs/mosaic-gui.log"
EOF
    echo -e "${GREEN}config.yaml created${NC}"
fi

# Final checks
echo ""
echo "================================================"
echo "   Installation Complete!"
echo "================================================"
echo ""
echo "To start the application:"
echo "  1. Activate the virtual environment:"
echo "     source venv/bin/activate"
echo ""
echo "  2. Run the Streamlit app:"
echo "     streamlit run app.py"
echo ""
echo "  3. Open your browser to:"
echo "     http://localhost:8501"
echo ""
echo "For command-line usage:"
echo "  python cli.py --smiles 'CCO' --format text"
echo "  python cli.py --name 'Aspirin' --output results.json"
echo ""
echo "For full MOSAIC model integration:"
echo "  1. Clone https://github.com/haoteli/MOSAIC"
echo "  2. Follow their setup instructions"
echo "  3. Copy model weights to ./Model_Weights/"
echo ""
echo -e "${GREEN}Setup complete!${NC}"
