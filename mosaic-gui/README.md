# MOSAIC-GUI: AI-Driven Retrosynthesis Planning Interface

<div align="center">

![MOSAIC](https://img.shields.io/badge/MOSAIC-Retrosynthesis-blue?style=for-the-badge)
![Python](https://img.shields.io/badge/Python-3.11+-green?style=for-the-badge)
![Streamlit](https://img.shields.io/badge/Streamlit-1.28+-red?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)

**A user-friendly GUI for AI-driven synthetic route planning with literature citations**

[Features](#features) • [Installation](#installation) • [Usage](#usage) • [Integration](#mosaic-integration) • [Citation](#citation)

</div>

---

## 🎯 Overview

MOSAIC-GUI provides an intuitive graphical interface for the [MOSAIC framework](https://github.com/haoteli/MOSAIC) (Multiple Optimized Specialists for AI-Driven Chemical Predictions), enabling researchers to:

- **Input molecules** via SMILES notation, compound names (PubChem), or CID
- **Generate step-by-step** retrosynthetic routes
- **View cited literature** for each reaction transformation
- **Export reports** in multiple formats (PDF, JSON, BibTeX)

## ✨ Features

### Input Options
| Input Type | Description | Example |
|------------|-------------|---------|
| **SMILES** | Direct molecular structure notation | `CC(=O)Oc1ccccc1C(=O)O` |
| **PubChem Name** | Common or IUPAC names | `Aspirin`, `Ibuprofen` |
| **PubChem CID** | PubChem Compound ID | `2244` |

### Output Features
- 📊 **Visual Molecule Rendering** - High-quality 2D structure images
- 🔄 **Step-by-Step Routes** - Clear retrosynthetic pathway
- 📚 **Literature Citations** - DOI-linked references for each transformation
- 📈 **Confidence Scores** - AI-generated reliability metrics
- 📝 **Reaction Conditions** - Temperature, solvent, reagent suggestions
- 💾 **Export Options** - PDF reports, BibTeX citations, JSON data

## 🚀 Installation

### Quick Start (Demo Mode)

```bash
# Clone this repository
git clone https://github.com/your-username/mosaic-gui.git
cd mosaic-gui

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run the application
streamlit run app.py
```

### Full Installation (With MOSAIC Models)

For full retrosynthesis prediction capability, you need the trained MOSAIC models:

1. **Clone the MOSAIC repository:**
```bash
git clone https://github.com/haoteli/MOSAIC.git
```

2. **Download trained models:**
   - Follow instructions in the [MOSAIC repository](https://github.com/haoteli/MOSAIC)
   - Download the 2,498 expert model weights
   - Download the FAISS index for Voronoi region assignment
   - Download the Kernel Metric Network weights

3. **Configure paths:**
```python
# In config.yaml or environment variables
MODEL_BASE_PATH=/path/to/Model_Weights
FAISS_INDEX_PATH=/path/to/faiss_index
KERNEL_NETWORK_PATH=/path/to/kernel_metric_network.pt
```

4. **Install CUDA dependencies (for GPU acceleration):**
```bash
pip install torch --index-url https://download.pytorch.org/whl/cu118
pip install faiss-gpu
```

## 📖 Usage

### Basic Usage

1. **Start the application:**
```bash
streamlit run app.py
```

2. **Open your browser** to `http://localhost:8501`

3. **Enter your target molecule:**
   - Type SMILES directly
   - Or search by compound name
   - Or enter PubChem CID

4. **Click "Plan Synthesis"** to generate the retrosynthetic route

5. **Review results:**
   - Examine each step's reaction type and conditions
   - Check confidence scores
   - Access literature citations
   - Export as needed

### Command Line Interface

```bash
# Single molecule prediction
python cli.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --output report.json

# Batch processing
python cli.py --input molecules.csv --output results/

# With custom parameters
python cli.py --smiles "CCO" --max-steps 15 --beam-width 10
```

### Python API

```python
from mosaic_interface import MOSAICEngine, MOSAICConfig

# Initialize engine
config = MOSAICConfig(
    model_base_path="./Model_Weights",
    max_steps=10,
    beam_width=5
)
engine = MOSAICEngine(config)
engine.initialize()

# Predict route
route = engine.plan_synthesis("CC(=O)Oc1ccccc1C(=O)O")

# Access results
print(f"Total steps: {route.num_steps}")
print(f"Estimated yield: {route.estimated_yield:.1f}%")

for step in route.steps:
    print(f"Step {step.step_number}: {step.reaction_type}")
    print(f"  Reactants: {step.reactants}")
    print(f"  Confidence: {step.confidence:.2%}")
```

## 🔬 MOSAIC Integration

### Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                      MOSAIC-GUI                              │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐       │
│  │   Streamlit  │  │   PubChem    │  │   Citation   │       │
│  │   Interface  │  │   API        │  │   Manager    │       │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘       │
│         │                 │                 │               │
│  ┌──────┴─────────────────┴─────────────────┴───────┐       │
│  │              MOSAIC Interface Layer               │       │
│  └──────────────────────┬────────────────────────────┘       │
└─────────────────────────┼────────────────────────────────────┘
                          │
┌─────────────────────────┼────────────────────────────────────┐
│                   MOSAIC Core                                │
│  ┌──────────────┐  ┌────┴─────┐  ┌──────────────────┐       │
│  │   Kernel     │  │  Voronoi │  │   2,498 Expert   │       │
│  │   Metric     │  │  FAISS   │  │   LLM Models     │       │
│  │   Network    │  │  Index   │  │   (Llama 3.1)    │       │
│  └──────────────┘  └──────────┘  └──────────────────┘       │
└──────────────────────────────────────────────────────────────┘
```

### Voronoi Region Assignment

MOSAIC partitions chemical space into 2,498 Voronoi regions. Each region is served by a specialized expert model:

1. **Molecular Fingerprint** - Compute Morgan fingerprint (2048 bits)
2. **Kernel Transformation** - Apply learned metric network
3. **Region Assignment** - Query FAISS index for nearest centroid
4. **Expert Selection** - Load appropriate specialist model
5. **Prediction** - Generate retrosynthesis with fine-tuned LLM

### Extending with Custom Experts

```python
from mosaic_interface import ExpertModel, MOSAICEngine

class CustomExpert(ExpertModel):
    """Custom expert for specific reaction types."""
    
    def predict(self, product_smiles, **kwargs):
        # Your custom prediction logic
        pass

# Register custom expert
engine = MOSAICEngine()
engine.experts[custom_id] = CustomExpert(
    expert_id=custom_id,
    model_path="./custom_models"
)
```

## 📊 Example Output

### Aspirin Retrosynthesis

```
Target: CC(=O)Oc1ccccc1C(=O)O (Aspirin)

Step 1: Esterification
├── Reactants: Salicylic acid + Acetic anhydride
├── Conditions: Reflux, H2SO4 catalyst
├── Confidence: 94%
├── Est. Yield: ~85%
└── Citations: Fischer, E. Ber. 1895 (DOI: 10.1002/cber.18950280321)

Step 2: Carboxylation (Kolbe-Schmitt)
├── Reactants: Phenol + CO2
├── Conditions: NaOH, 125°C, 100 atm
├── Confidence: 88%
├── Est. Yield: ~70%
└── Citations: Kolbe, H. Ann. Chem. 1860 (DOI: 10.1002/jlac.18601130102)

Overall Yield Estimate: ~60%
Total Confidence: 83%
```

## 🎓 Citation

If you use MOSAIC-GUI in your research, please cite:

```bibtex
@article{mosaic2024,
  title={Collective Intelligence for AI-Assisted Chemical Synthesis},
  author={Li, Haote and others},
  journal={},
  year={2024},
  url={https://github.com/haoteli/MOSAIC}
}
```

## 📁 Project Structure

```
mosaic-gui/
├── app.py                  # Main Streamlit application
├── mosaic_interface.py     # MOSAIC model interface
├── requirements.txt        # Python dependencies
├── README.md              # This file
├── config.yaml            # Configuration (optional)
├── cli.py                 # Command-line interface
└── utils/
    ├── pubchem.py         # PubChem API utilities
    ├── citation.py        # Citation management
    └── export.py          # Report generation
```

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- [MOSAIC](https://github.com/haoteli/MOSAIC) - The underlying AI framework
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - Chemical compound database
- [RDKit](https://www.rdkit.org/) - Cheminformatics toolkit
- [Streamlit](https://streamlit.io/) - Web application framework

---

<div align="center">

**Made with ❤️ for the chemistry research community**

[Report Bug](https://github.com/your-username/mosaic-gui/issues) • [Request Feature](https://github.com/your-username/mosaic-gui/issues)

</div>
