---
title: MOSAIC Retrosynthesis Planner
emoji: 🧪
colorFrom: green
colorTo: blue
sdk: streamlit
sdk_version: 1.31.0
app_file: app.py
pinned: true
license: mit
tags:
  - chemistry
  - retrosynthesis
  - drug-discovery
  - ai
  - llm
  - natural-products
short_description: AI-powered retrosynthetic route planning with literature citations
---

# MOSAIC-GUI: AI-Powered Retrosynthesis Planner

🧪 **MOSAIC** (Multiple Optimized Specialists for AI-Driven Chemical predictions) is a computational framework that provides step-by-step synthetic routes for target molecules.

## Features

- 🔬 **SMILES Input**: Enter molecules using SMILES notation
- 🔍 **PubChem Lookup**: Search by compound name or CID
- 📊 **Step-by-Step Routes**: Detailed synthetic pathway with conditions
- 📚 **Literature Citations**: DOI-linked references for each transformation
- 📈 **Confidence Scores**: AI-predicted success probability
- 📥 **Export Options**: PDF, JSON, CSV, BibTeX

## Usage

1. Select input mode (SMILES, Name, or CID)
2. Enter your target molecule
3. Click "Generate Retrosynthetic Route"
4. Review the step-by-step synthesis plan
5. Export results in your preferred format

## Example Molecules

- **Aspirin**: `CC(=O)Oc1ccccc1C(=O)O`
- **Ibuprofen**: `CC(C)Cc1ccc(C(C)C(=O)O)cc1`
- **Caffeine**: `Cn1cnc2c1c(=O)n(c(=O)n2C)C`

## Citation

If you use this tool in your research, please cite:

```bibtex
@article{mosaic2024,
  title={Collective Intelligence for AI-Assisted Chemical Synthesis},
  author={Li, Haote and others},
  journal={TBD},
  year={2024}
}
```

## License

MIT License - Free for academic and commercial use.

## Contact

For questions or collaborations, please open an issue on GitHub.
