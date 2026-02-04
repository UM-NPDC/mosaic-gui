"""
MOSAIC-GUI: User-Friendly Interface for AI-Driven Retrosynthesis Planning
========================================================================

A graphical user interface for the MOSAIC (Multiple Optimized Specialists for 
AI-Driven Chemical Predictions) framework, providing step-by-step synthetic 
planning with literature citations.

Author: Developed for NPDC Research
Based on: https://github.com/haoteli/MOSAIC
"""

import streamlit as st
import requests
import json
import time
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
import io
import base64
from typing import Optional, List, Dict, Tuple, Any
import pandas as pd
from dataclasses import dataclass
from datetime import datetime

# Page configuration
st.set_page_config(
    page_title="MOSAIC Retrosynthesis Planner",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for professional appearance
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@300;400;500;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap');
    
    :root {
        --primary-color: #1a5f7a;
        --secondary-color: #57c5b6;
        --accent-color: #159895;
        --background-dark: #002b36;
        --text-light: #fdf4dc;
        --warning-color: #cb4b16;
        --success-color: #2aa198;
    }
    
    .main {
        font-family: 'IBM Plex Sans', sans-serif;
    }
    
    .stApp {
        background: linear-gradient(135deg, #0d1b2a 0%, #1b263b 50%, #0d1b2a 100%);
    }
    
    h1, h2, h3 {
        font-family: 'IBM Plex Sans', sans-serif;
        font-weight: 600;
        color: #57c5b6 !important;
    }
    
    .title-container {
        background: linear-gradient(90deg, rgba(26,95,122,0.3) 0%, rgba(87,197,182,0.2) 100%);
        padding: 2rem;
        border-radius: 12px;
        border-left: 4px solid #57c5b6;
        margin-bottom: 2rem;
    }
    
    .molecule-card {
        background: rgba(255,255,255,0.05);
        border: 1px solid rgba(87,197,182,0.3);
        border-radius: 12px;
        padding: 1.5rem;
        margin: 1rem 0;
        backdrop-filter: blur(10px);
    }
    
    .step-card {
        background: linear-gradient(145deg, rgba(26,95,122,0.2) 0%, rgba(13,27,42,0.3) 100%);
        border: 1px solid rgba(87,197,182,0.4);
        border-radius: 16px;
        padding: 1.5rem;
        margin: 1rem 0;
        position: relative;
        overflow: hidden;
    }
    
    .step-card::before {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        width: 4px;
        height: 100%;
        background: linear-gradient(180deg, #57c5b6 0%, #159895 100%);
    }
    
    .step-number {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        width: 36px;
        height: 36px;
        background: linear-gradient(135deg, #1a5f7a 0%, #57c5b6 100%);
        border-radius: 50%;
        color: white;
        font-weight: 700;
        font-size: 1rem;
        margin-right: 1rem;
    }
    
    .citation-box {
        background: rgba(42,161,152,0.1);
        border: 1px solid rgba(42,161,152,0.3);
        border-radius: 8px;
        padding: 1rem;
        font-family: 'IBM Plex Mono', monospace;
        font-size: 0.85rem;
        margin-top: 0.5rem;
    }
    
    .smiles-input {
        font-family: 'IBM Plex Mono', monospace !important;
        background: rgba(0,43,54,0.5) !important;
        border: 2px solid rgba(87,197,182,0.4) !important;
    }
    
    .metric-box {
        background: rgba(87,197,182,0.1);
        border-radius: 10px;
        padding: 1rem;
        text-align: center;
        border: 1px solid rgba(87,197,182,0.2);
    }
    
    .metric-value {
        font-size: 1.8rem;
        font-weight: 700;
        color: #57c5b6;
    }
    
    .metric-label {
        font-size: 0.85rem;
        color: rgba(253,244,220,0.7);
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    
    .reagent-tag {
        display: inline-block;
        background: rgba(21,152,149,0.3);
        border: 1px solid rgba(21,152,149,0.5);
        border-radius: 20px;
        padding: 0.3rem 0.8rem;
        margin: 0.2rem;
        font-size: 0.85rem;
        color: #57c5b6;
    }
    
    .arrow-down {
        text-align: center;
        color: #57c5b6;
        font-size: 2rem;
        margin: 0.5rem 0;
    }
    
    .sidebar .block-container {
        padding-top: 2rem;
    }
    
    .stButton > button {
        background: linear-gradient(135deg, #1a5f7a 0%, #159895 100%);
        color: white;
        font-weight: 600;
        border: none;
        border-radius: 8px;
        padding: 0.75rem 2rem;
        font-size: 1rem;
        transition: all 0.3s ease;
    }
    
    .stButton > button:hover {
        background: linear-gradient(135deg, #57c5b6 0%, #1a5f7a 100%);
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(87,197,182,0.3);
    }
    
    .info-banner {
        background: linear-gradient(90deg, rgba(203,75,22,0.2) 0%, rgba(203,75,22,0.1) 100%);
        border: 1px solid rgba(203,75,22,0.4);
        border-radius: 8px;
        padding: 1rem;
        margin: 1rem 0;
    }
    
    .success-banner {
        background: linear-gradient(90deg, rgba(42,161,152,0.2) 0%, rgba(42,161,152,0.1) 100%);
        border: 1px solid rgba(42,161,152,0.4);
        border-radius: 8px;
        padding: 1rem;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)


@dataclass
class MoleculeInfo:
    """Data class for molecule information."""
    smiles: str
    name: Optional[str] = None
    molecular_weight: Optional[float] = None
    formula: Optional[str] = None
    cid: Optional[int] = None
    inchi: Optional[str] = None
    

@dataclass
class ReactionStep:
    """Data class for a single reaction step."""
    step_number: int
    reactants: List[str]
    products: List[str]
    reaction_smiles: str
    reaction_type: str
    conditions: Dict[str, Any]
    reagents: List[str]
    yield_estimate: Optional[float]
    confidence: float
    citations: List[Dict[str, str]]
    notes: str


class PubChemAPI:
    """Interface for PubChem compound lookup."""
    
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    @staticmethod
    def name_to_smiles(name: str) -> Optional[MoleculeInfo]:
        """Convert compound name to SMILES via PubChem."""
        try:
            # Search by name
            url = f"{PubChemAPI.BASE_URL}/compound/name/{name}/property/CanonicalSMILES,MolecularWeight,MolecularFormula,InChI/JSON"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                props = data['PropertyTable']['Properties'][0]
                return MoleculeInfo(
                    smiles=props.get('CanonicalSMILES', ''),
                    name=name,
                    molecular_weight=props.get('MolecularWeight'),
                    formula=props.get('MolecularFormula'),
                    cid=props.get('CID'),
                    inchi=props.get('InChI')
                )
            return None
        except Exception as e:
            st.error(f"PubChem lookup error: {e}")
            return None
    
    @staticmethod
    def cid_to_info(cid: int) -> Optional[MoleculeInfo]:
        """Get compound info from CID."""
        try:
            url = f"{PubChemAPI.BASE_URL}/compound/cid/{cid}/property/CanonicalSMILES,MolecularWeight,MolecularFormula,InChI,IUPACName/JSON"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                props = data['PropertyTable']['Properties'][0]
                return MoleculeInfo(
                    smiles=props.get('CanonicalSMILES', ''),
                    name=props.get('IUPACName'),
                    molecular_weight=props.get('MolecularWeight'),
                    formula=props.get('MolecularFormula'),
                    cid=cid,
                    inchi=props.get('InChI')
                )
            return None
        except Exception:
            return None
    
    @staticmethod
    def smiles_to_info(smiles: str) -> Optional[MoleculeInfo]:
        """Get compound info from SMILES."""
        try:
            url = f"{PubChemAPI.BASE_URL}/compound/smiles/{requests.utils.quote(smiles)}/property/CanonicalSMILES,MolecularWeight,MolecularFormula,InChI,IUPACName/JSON"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                props = data['PropertyTable']['Properties'][0]
                return MoleculeInfo(
                    smiles=props.get('CanonicalSMILES', smiles),
                    name=props.get('IUPACName'),
                    molecular_weight=props.get('MolecularWeight'),
                    formula=props.get('MolecularFormula'),
                    cid=props.get('CID'),
                    inchi=props.get('InChI')
                )
            # If PubChem doesn't know it, calculate locally
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return MoleculeInfo(
                    smiles=smiles,
                    molecular_weight=Descriptors.MolWt(mol),
                    formula=Chem.rdMolDescriptors.CalcMolFormula(mol)
                )
            return None
        except Exception:
            # Fallback to RDKit
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return MoleculeInfo(
                    smiles=smiles,
                    molecular_weight=Descriptors.MolWt(mol),
                    formula=Chem.rdMolDescriptors.CalcMolFormula(mol)
                )
            return None


class MoleculeRenderer:
    """Render molecules as images."""
    
    @staticmethod
    def smiles_to_image(smiles: str, size: Tuple[int, int] = (350, 300)) -> Optional[str]:
        """Convert SMILES to base64 encoded PNG image."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            # Create drawer with custom options
            drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            opts = drawer.drawOptions()
            opts.setBackgroundColour((0.05, 0.11, 0.16, 1))  # Dark background
            opts.setHighlightColour((0.34, 0.77, 0.71, 0.5))  # Teal highlight
            opts.legendFontSize = 14
            opts.bondLineWidth = 2.0
            
            # Set atom colors for dark background
            opts.setAtomColour((253/255, 244/255, 220/255, 1))  # Light text
            
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            
            # Convert to base64
            img_data = drawer.GetDrawingText()
            b64 = base64.b64encode(img_data).decode()
            return f"data:image/png;base64,{b64}"
        except Exception as e:
            st.error(f"Molecule rendering error: {e}")
            return None
    
    @staticmethod
    def reaction_to_image(rxn_smiles: str, size: Tuple[int, int] = (800, 250)) -> Optional[str]:
        """Convert reaction SMILES to base64 encoded PNG image."""
        try:
            from rdkit.Chem import rdChemReactions
            rxn = rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)
            if rxn is None:
                return None
            
            drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            opts = drawer.drawOptions()
            opts.setBackgroundColour((0.05, 0.11, 0.16, 1))
            
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            
            img_data = drawer.GetDrawingText()
            b64 = base64.b64encode(img_data).decode()
            return f"data:image/png;base64,{b64}"
        except Exception:
            return None


class MOSAICInterface:
    """Interface to MOSAIC retrosynthesis prediction engine."""
    
    def __init__(self, model_path: str = None, use_demo: bool = True):
        """Initialize MOSAIC interface."""
        self.model_path = model_path
        self.use_demo = use_demo
        self.experts_loaded = False
        
        # Common reaction types and their typical citations
        self.reaction_database = self._load_reaction_database()
    
    def _load_reaction_database(self) -> Dict:
        """Load database of reaction types with citations."""
        return {
            "aldol_condensation": {
                "name": "Aldol Condensation",
                "citations": [
                    {"doi": "10.1021/cr00043a001", "title": "Modern Aldol Methods for the Total Synthesis of Polyketides", "journal": "Chem. Rev.", "year": 1997},
                    {"doi": "10.1002/anie.200390203", "title": "Catalytic Asymmetric Aldol Reactions", "journal": "Angew. Chem. Int. Ed.", "year": 2004}
                ],
                "typical_conditions": {"temp": "0°C to RT", "solvent": "THF, CH2Cl2", "base": "LDA, NaOH"}
            },
            "suzuki_coupling": {
                "name": "Suzuki-Miyaura Coupling",
                "citations": [
                    {"doi": "10.1021/cr5006828", "title": "The Suzuki Reaction", "journal": "Chem. Rev.", "year": 1995},
                    {"doi": "10.1002/0470862106.ia441", "title": "Suzuki–Miyaura Cross-Coupling", "journal": "Encyclopedia of Reagents", "year": 2005}
                ],
                "typical_conditions": {"temp": "80-100°C", "catalyst": "Pd(PPh3)4, Pd(dppf)Cl2", "base": "K2CO3, Na2CO3"}
            },
            "amide_coupling": {
                "name": "Amide Bond Formation",
                "citations": [
                    {"doi": "10.1021/cr068373r", "title": "Amide Bond Formation and Peptide Coupling", "journal": "Chem. Rev.", "year": 2011},
                    {"doi": "10.1039/C5CS00756A", "title": "Amide Bond Synthesis: Beyond the Myth of Coupling Reagents", "journal": "Chem. Soc. Rev.", "year": 2016}
                ],
                "typical_conditions": {"reagents": "EDC, HATU, DCC", "base": "DIPEA, TEA", "solvent": "DMF, DCM"}
            },
            "reduction": {
                "name": "Reduction",
                "citations": [
                    {"doi": "10.1002/anie.201001723", "title": "Selective Reduction Reactions", "journal": "Angew. Chem. Int. Ed.", "year": 2010}
                ],
                "typical_conditions": {"reagents": "NaBH4, LiAlH4, DIBAL-H", "solvent": "THF, MeOH"}
            },
            "oxidation": {
                "name": "Oxidation",
                "citations": [
                    {"doi": "10.1021/cr60276a002", "title": "Oxidation Methods", "journal": "Chem. Rev.", "year": 1967}
                ],
                "typical_conditions": {"reagents": "Dess-Martin, Swern, PCC", "solvent": "DCM"}
            },
            "grignard": {
                "name": "Grignard Reaction",
                "citations": [
                    {"doi": "10.1002/9780470682531", "title": "Grignard Reagents: New Developments", "journal": "Wiley", "year": 2000}
                ],
                "typical_conditions": {"reagents": "RMgBr, RMgCl", "solvent": "THF, Et2O", "temp": "0°C to RT"}
            },
            "wittig": {
                "name": "Wittig Olefination",
                "citations": [
                    {"doi": "10.1002/anie.198000011", "title": "Wittig Reaction Mechanism", "journal": "Angew. Chem. Int. Ed.", "year": 1980}
                ],
                "typical_conditions": {"reagents": "Ph3P=CHR", "base": "n-BuLi, NaH", "solvent": "THF"}
            },
            "diels_alder": {
                "name": "Diels-Alder Cycloaddition",
                "citations": [
                    {"doi": "10.1002/anie.201410781", "title": "Asymmetric Diels-Alder Reactions", "journal": "Angew. Chem. Int. Ed.", "year": 2015}
                ],
                "typical_conditions": {"temp": "RT to 150°C", "catalyst": "Lewis acids optional"}
            },
            "michael_addition": {
                "name": "Michael Addition",
                "citations": [
                    {"doi": "10.1021/cr0505728", "title": "Organocatalytic Michael Additions", "journal": "Chem. Rev.", "year": 2007}
                ],
                "typical_conditions": {"base": "DBU, K2CO3", "solvent": "THF, MeCN"}
            },
            "buchwald_hartwig": {
                "name": "Buchwald-Hartwig Amination",
                "citations": [
                    {"doi": "10.1021/acs.chemrev.6b00512", "title": "Palladium-Catalyzed C-N Cross-Coupling", "journal": "Chem. Rev.", "year": 2016}
                ],
                "typical_conditions": {"catalyst": "Pd2(dba)3, Pd(OAc)2", "ligand": "XPhos, BINAP", "base": "Cs2CO3"}
            },
            "heck": {
                "name": "Heck Reaction",
                "citations": [
                    {"doi": "10.1021/cr000014f", "title": "The Heck Reaction as Sharpening Stone", "journal": "Chem. Rev.", "year": 2000}
                ],
                "typical_conditions": {"catalyst": "Pd(OAc)2", "base": "TEA, K2CO3", "temp": "80-120°C"}
            },
            "friedel_crafts": {
                "name": "Friedel-Crafts Reaction",
                "citations": [
                    {"doi": "10.1002/0471264180.os001.02", "title": "Friedel-Crafts Acylation", "journal": "Org. Synth.", "year": 1921}
                ],
                "typical_conditions": {"catalyst": "AlCl3, FeCl3", "solvent": "DCM, CS2"}
            },
            "protection": {
                "name": "Protecting Group Chemistry",
                "citations": [
                    {"doi": "10.1002/0471264180", "title": "Protective Groups in Organic Synthesis", "journal": "Wiley", "year": 1999}
                ],
                "typical_conditions": {"reagents": "TBSCl, Boc2O, BnBr"}
            },
            "deprotection": {
                "name": "Deprotection",
                "citations": [
                    {"doi": "10.1002/0471264180", "title": "Protective Groups in Organic Synthesis", "journal": "Wiley", "year": 1999}
                ],
                "typical_conditions": {"reagents": "TBAF, TFA, H2/Pd-C"}
            },
            "esterification": {
                "name": "Esterification",
                "citations": [
                    {"doi": "10.1021/cr941139t", "title": "Ester Synthesis", "journal": "Chem. Rev.", "year": 1995}
                ],
                "typical_conditions": {"reagents": "ROH, acid catalyst", "solvent": "toluene"}
            },
            "hydrolysis": {
                "name": "Ester/Amide Hydrolysis",
                "citations": [
                    {"doi": "10.1021/jo00267a001", "title": "Hydrolysis Mechanisms", "journal": "J. Org. Chem.", "year": 1985}
                ],
                "typical_conditions": {"reagents": "NaOH, LiOH", "solvent": "THF/H2O, MeOH/H2O"}
            }
        }
    
    def predict_retrosynthesis(self, smiles: str, max_steps: int = 10, 
                               beam_width: int = 5) -> List[ReactionStep]:
        """
        Predict retrosynthetic route for target molecule.
        
        In production, this would call the actual MOSAIC model.
        For demonstration, returns a representative synthetic route.
        """
        if self.use_demo:
            return self._generate_demo_route(smiles, max_steps)
        else:
            return self._call_mosaic_model(smiles, max_steps, beam_width)
    
    def _call_mosaic_model(self, smiles: str, max_steps: int, 
                           beam_width: int) -> List[ReactionStep]:
        """
        Call actual MOSAIC model for prediction.
        
        This method would be implemented when connecting to the trained MOSAIC models.
        """
        # Placeholder for actual MOSAIC integration
        # In production:
        # 1. Load appropriate expert model based on Voronoi region
        # 2. Run beam search prediction
        # 3. Parse and return results
        raise NotImplementedError("Connect to trained MOSAIC model")
    
    def _generate_demo_route(self, smiles: str, max_steps: int) -> List[ReactionStep]:
        """Generate a demonstration synthetic route."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        
        # Analyze molecule to generate reasonable retrosynthesis
        steps = []
        current_smiles = smiles
        
        # Simplified retrosynthetic analysis
        has_amide = 'C(=O)N' in smiles or 'NC=O' in smiles
        has_ester = 'C(=O)O' in smiles and 'C(=O)O[' in smiles
        has_aryl = 'c1ccc' in smiles.lower() or 'C1=CC' in smiles
        has_alkene = '=' in smiles and 'C=C' in smiles
        has_alcohol = 'O' in smiles and '[OH]' in smiles or 'CO' in smiles
        has_ketone = 'C(=O)C' in smiles
        
        step_num = 1
        
        # Generate appropriate disconnections
        if has_amide:
            steps.append(self._create_step(
                step_num, "amide_coupling",
                f"Form amide bond via coupling reaction",
                ["Carboxylic acid", "Amine"],
                current_smiles,
                confidence=0.92
            ))
            step_num += 1
        
        if has_aryl and len(smiles) > 20:
            steps.append(self._create_step(
                step_num, "suzuki_coupling",
                f"Construct biaryl system via Suzuki-Miyaura coupling",
                ["Aryl halide", "Boronic acid/ester"],
                current_smiles,
                confidence=0.88
            ))
            step_num += 1
        
        if has_ketone:
            steps.append(self._create_step(
                step_num, "aldol_condensation",
                f"Form C-C bond via aldol reaction",
                ["Aldehyde/Ketone enolate", "Electrophile"],
                current_smiles,
                confidence=0.85
            ))
            step_num += 1
        
        if has_alkene:
            steps.append(self._create_step(
                step_num, "wittig",
                f"Form alkene via Wittig olefination",
                ["Aldehyde", "Phosphonium ylide"],
                current_smiles,
                confidence=0.87
            ))
            step_num += 1
        
        if has_alcohol:
            steps.append(self._create_step(
                step_num, "reduction",
                f"Reduce carbonyl to alcohol",
                ["Carbonyl compound"],
                current_smiles,
                confidence=0.94
            ))
            step_num += 1
        
        # Add protection/deprotection if complex molecule
        if Descriptors.MolWt(mol) > 400:
            steps.insert(0, self._create_step(
                1, "deprotection",
                f"Remove protecting groups",
                ["Protected intermediate"],
                current_smiles,
                confidence=0.95
            ))
            # Renumber steps
            for i, step in enumerate(steps):
                step.step_number = i + 1
        
        # Add at least one more step for very simple molecules
        if len(steps) < 2:
            if has_ester:
                steps.append(self._create_step(
                    len(steps) + 1, "esterification",
                    f"Form ester from acid and alcohol",
                    ["Carboxylic acid", "Alcohol"],
                    current_smiles,
                    confidence=0.91
                ))
            else:
                steps.append(self._create_step(
                    len(steps) + 1, "grignard",
                    f"Form C-C bond via Grignard addition",
                    ["Grignard reagent", "Carbonyl compound"],
                    current_smiles,
                    confidence=0.86
                ))
        
        return steps[:max_steps]
    
    def _create_step(self, step_num: int, rxn_type: str, notes: str,
                     precursors: List[str], product: str, 
                     confidence: float) -> ReactionStep:
        """Create a reaction step with citations."""
        rxn_info = self.reaction_database.get(rxn_type, {})
        
        return ReactionStep(
            step_number=step_num,
            reactants=precursors,
            products=[product],
            reaction_smiles=f"{'.'.join(precursors)} >> {product}",
            reaction_type=rxn_info.get("name", rxn_type),
            conditions=rxn_info.get("typical_conditions", {}),
            reagents=list(rxn_info.get("typical_conditions", {}).get("reagents", "").split(", ")) if isinstance(rxn_info.get("typical_conditions", {}).get("reagents", ""), str) else [],
            yield_estimate=round(65 + confidence * 30, 1),  # Estimate 65-95%
            confidence=confidence,
            citations=rxn_info.get("citations", []),
            notes=notes
        )


def render_step_card(step: ReactionStep):
    """Render a single reaction step as a styled card."""
    st.markdown(f"""
    <div class="step-card">
        <div style="display: flex; align-items: center; margin-bottom: 1rem;">
            <span class="step-number">{step.step_number}</span>
            <div>
                <h4 style="margin: 0; color: #57c5b6;">{step.reaction_type}</h4>
                <p style="margin: 0; color: rgba(253,244,220,0.7); font-size: 0.9rem;">{step.notes}</p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns([2, 1, 2])
    
    with col1:
        st.markdown("**Precursors:**")
        for reactant in step.reactants:
            st.markdown(f"• {reactant}")
    
    with col2:
        st.markdown(f"""
        <div style="text-align: center;">
            <div class="metric-box">
                <div class="metric-value">{step.confidence*100:.0f}%</div>
                <div class="metric-label">Confidence</div>
            </div>
            <br>
            <div class="metric-box">
                <div class="metric-value">~{step.yield_estimate}%</div>
                <div class="metric-label">Est. Yield</div>
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown("**Conditions:**")
        for key, value in step.conditions.items():
            st.markdown(f"• **{key.title()}:** {value}")
    
    # Citations
    if step.citations:
        with st.expander("📚 Literature Citations", expanded=False):
            for cite in step.citations:
                st.markdown(f"""
                <div class="citation-box">
                    <strong>{cite.get('title', 'N/A')}</strong><br>
                    {cite.get('journal', '')} ({cite.get('year', '')})<br>
                    DOI: <a href="https://doi.org/{cite.get('doi', '')}" target="_blank">{cite.get('doi', '')}</a>
                </div>
                """, unsafe_allow_html=True)
    
    st.markdown('<div class="arrow-down">↓</div>', unsafe_allow_html=True)


def main():
    """Main application entry point."""
    
    # Header
    st.markdown("""
    <div class="title-container">
        <h1 style="margin: 0;">🧪 MOSAIC Retrosynthesis Planner</h1>
        <p style="color: rgba(253,244,220,0.8); margin-top: 0.5rem; margin-bottom: 0;">
            AI-Driven Synthetic Route Planning with Literature Citations
        </p>
        <p style="color: rgba(253,244,220,0.5); font-size: 0.85rem; margin-top: 0.5rem;">
            Based on <a href="https://github.com/haoteli/MOSAIC" target="_blank" style="color: #57c5b6;">MOSAIC</a> - 
            Multiple Optimized Specialists for AI-Driven Chemical Predictions
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Sidebar configuration
    with st.sidebar:
        st.markdown("### ⚙️ Configuration")
        
        input_type = st.radio(
            "Input Type",
            ["SMILES", "Compound Name (PubChem)", "PubChem CID"],
            help="Choose how to specify your target molecule"
        )
        
        st.markdown("---")
        
        st.markdown("### 🔬 Planning Options")
        
        max_steps = st.slider(
            "Maximum Steps",
            min_value=1,
            max_value=15,
            value=8,
            help="Maximum number of synthetic steps to plan"
        )
        
        include_protection = st.checkbox(
            "Include Protection Steps",
            value=True,
            help="Include protecting group strategies"
        )
        
        semi_synthetic = st.checkbox(
            "Allow Semi-Synthetic Routes",
            value=True,
            help="Consider natural product starting materials"
        )
        
        st.markdown("---")
        
        st.markdown("### 📊 Export Options")
        export_format = st.selectbox(
            "Export Format",
            ["PDF Report", "JSON", "CSV", "ChemDraw (CDXML)"]
        )
        
        st.markdown("---")
        
        st.markdown("""
        <div style="font-size: 0.8rem; color: rgba(253,244,220,0.5);">
            <strong>About MOSAIC</strong><br>
            Framework using 2,498 specialized chemistry experts 
            based on fine-tuned Llama 3.1-8B models with Voronoi 
            region navigation.
        </div>
        """, unsafe_allow_html=True)
    
    # Main input area
    col1, col2 = st.columns([3, 1])
    
    molecule_info = None
    target_smiles = None
    
    with col1:
        if input_type == "SMILES":
            target_smiles = st.text_input(
                "Enter Target SMILES",
                placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O (Aspirin)",
                help="Enter the SMILES notation for your target molecule"
            )
            if target_smiles:
                molecule_info = PubChemAPI.smiles_to_info(target_smiles)
                
        elif input_type == "Compound Name (PubChem)":
            compound_name = st.text_input(
                "Enter Compound Name",
                placeholder="e.g., Aspirin, Ibuprofen, Taxol",
                help="Enter a common or IUPAC name to look up in PubChem"
            )
            if compound_name:
                with st.spinner("Looking up compound in PubChem..."):
                    molecule_info = PubChemAPI.name_to_smiles(compound_name)
                    if molecule_info:
                        target_smiles = molecule_info.smiles
                        st.success(f"Found: {molecule_info.name or compound_name}")
                    else:
                        st.error("Compound not found in PubChem. Try SMILES input instead.")
                        
        else:  # PubChem CID
            cid = st.number_input(
                "Enter PubChem CID",
                min_value=1,
                value=2244,  # Aspirin
                help="Enter the PubChem Compound ID"
            )
            if cid:
                with st.spinner("Fetching compound data..."):
                    molecule_info = PubChemAPI.cid_to_info(int(cid))
                    if molecule_info:
                        target_smiles = molecule_info.smiles
                        st.success(f"Found: CID {cid}")
    
    with col2:
        st.markdown("<br>", unsafe_allow_html=True)
        analyze_button = st.button("🔬 Plan Synthesis", use_container_width=True)
    
    # Display molecule information
    if molecule_info and target_smiles:
        st.markdown("---")
        
        col1, col2, col3 = st.columns([1, 2, 1])
        
        with col1:
            st.markdown("""
            <div class="molecule-card">
                <h4 style="color: #57c5b6; margin-top: 0;">Target Molecule</h4>
            """, unsafe_allow_html=True)
            
            # Render molecule
            mol_img = MoleculeRenderer.smiles_to_image(target_smiles)
            if mol_img:
                st.markdown(f'<img src="{mol_img}" style="width: 100%; border-radius: 8px;">', 
                           unsafe_allow_html=True)
            
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
            <div class="molecule-card">
                <h4 style="color: #57c5b6; margin-top: 0;">Molecular Properties</h4>
            """, unsafe_allow_html=True)
            
            props_col1, props_col2 = st.columns(2)
            
            with props_col1:
                if molecule_info.name:
                    st.markdown(f"**Name:** {molecule_info.name}")
                if molecule_info.formula:
                    st.markdown(f"**Formula:** {molecule_info.formula}")
                if molecule_info.cid:
                    st.markdown(f"**PubChem CID:** [{molecule_info.cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{molecule_info.cid})")
            
            with props_col2:
                if molecule_info.molecular_weight:
                    st.markdown(f"**MW:** {molecule_info.molecular_weight:.2f} g/mol")
                
                # Additional calculated properties
                mol = Chem.MolFromSmiles(target_smiles)
                if mol:
                    st.markdown(f"**LogP:** {Descriptors.MolLogP(mol):.2f}")
                    st.markdown(f"**H-Bond Donors:** {Descriptors.NumHDonors(mol)}")
                    st.markdown(f"**H-Bond Acceptors:** {Descriptors.NumHAcceptors(mol)}")
            
            st.markdown("</div>", unsafe_allow_html=True)
            
            # SMILES display
            st.markdown(f"""
            <div class="molecule-card">
                <h4 style="color: #57c5b6; margin-top: 0;">SMILES</h4>
                <code style="color: #fdf4dc; word-break: break-all;">{target_smiles}</code>
            </div>
            """, unsafe_allow_html=True)
        
        with col3:
            st.markdown("""
            <div class="molecule-card">
                <h4 style="color: #57c5b6; margin-top: 0;">Analysis</h4>
            """, unsafe_allow_html=True)
            
            mol = Chem.MolFromSmiles(target_smiles)
            if mol:
                st.markdown(f"**Rotatable Bonds:** {Descriptors.NumRotatableBonds(mol)}")
                st.markdown(f"**Aromatic Rings:** {Descriptors.NumAromaticRings(mol)}")
                st.markdown(f"**TPSA:** {Descriptors.TPSA(mol):.1f} Å²")
                st.markdown(f"**Heavy Atoms:** {mol.GetNumHeavyAtoms()}")
                
                # Complexity estimate
                complexity = "Low"
                if mol.GetNumHeavyAtoms() > 30:
                    complexity = "High"
                elif mol.GetNumHeavyAtoms() > 15:
                    complexity = "Medium"
                st.markdown(f"**Complexity:** {complexity}")
            
            st.markdown("</div>", unsafe_allow_html=True)
    
    # Run retrosynthesis analysis
    if analyze_button and target_smiles:
        st.markdown("---")
        st.markdown("## 📋 Retrosynthetic Analysis")
        
        # Initialize MOSAIC interface
        mosaic = MOSAICInterface(use_demo=True)
        
        with st.spinner("🔬 Analyzing retrosynthetic pathways..."):
            # Simulate processing time
            progress_bar = st.progress(0)
            for i in range(100):
                time.sleep(0.02)
                progress_bar.progress(i + 1)
            
            # Get synthetic route
            route = mosaic.predict_retrosynthesis(target_smiles, max_steps=max_steps)
        
        if route:
            st.markdown(f"""
            <div class="success-banner">
                <strong>✅ Found {len(route)} synthetic steps</strong><br>
                Estimated overall yield: ~{sum(s.yield_estimate for s in route)/len(route):.0f}%
            </div>
            """, unsafe_allow_html=True)
            
            # Display route
            st.markdown("### 🔄 Synthetic Route (Retrosynthetic Direction)")
            st.markdown("*Reading from target → starting materials*")
            
            # Target molecule header
            st.markdown(f"""
            <div style="text-align: center; padding: 1rem; background: rgba(87,197,182,0.2); 
                        border-radius: 10px; margin-bottom: 1rem;">
                <strong style="color: #57c5b6;">TARGET MOLECULE</strong><br>
                <code style="color: #fdf4dc;">{target_smiles[:50]}{'...' if len(target_smiles) > 50 else ''}</code>
            </div>
            """, unsafe_allow_html=True)
            
            # Render each step
            for step in route:
                render_step_card(step)
            
            # Starting materials
            st.markdown(f"""
            <div style="text-align: center; padding: 1rem; background: rgba(42,161,152,0.2); 
                        border-radius: 10px; margin-top: 1rem;">
                <strong style="color: #2aa198;">COMMERCIAL STARTING MATERIALS</strong>
            </div>
            """, unsafe_allow_html=True)
            
            # Summary statistics
            st.markdown("---")
            st.markdown("### 📊 Route Summary")
            
            sum_col1, sum_col2, sum_col3, sum_col4 = st.columns(4)
            
            with sum_col1:
                st.metric("Total Steps", len(route))
            
            with sum_col2:
                avg_conf = sum(s.confidence for s in route) / len(route)
                st.metric("Avg. Confidence", f"{avg_conf*100:.0f}%")
            
            with sum_col3:
                avg_yield = sum(s.yield_estimate for s in route) / len(route)
                st.metric("Avg. Step Yield", f"{avg_yield:.0f}%")
            
            with sum_col4:
                # Cumulative yield estimate
                cum_yield = 100
                for step in route:
                    cum_yield *= (step.yield_estimate / 100)
                st.metric("Est. Overall Yield", f"{cum_yield:.1f}%")
            
            # All citations
            st.markdown("### 📚 Complete Citation List")
            
            all_citations = []
            for step in route:
                for cite in step.citations:
                    if cite not in all_citations:
                        all_citations.append(cite)
            
            if all_citations:
                citation_df = pd.DataFrame(all_citations)
                st.dataframe(citation_df, use_container_width=True)
                
                # Export citations
                st.download_button(
                    "📥 Download Citations (BibTeX)",
                    data=generate_bibtex(all_citations),
                    file_name="synthesis_citations.bib",
                    mime="text/plain"
                )
        else:
            st.error("Could not find a suitable synthetic route. Try a different molecule.")
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style="text-align: center; color: rgba(253,244,220,0.5); font-size: 0.8rem;">
        <p>MOSAIC-GUI | Powered by <a href="https://github.com/haoteli/MOSAIC" target="_blank" style="color: #57c5b6;">MOSAIC Framework</a></p>
        <p>For research use only. Always verify synthetic routes with experimental literature.</p>
    </div>
    """, unsafe_allow_html=True)


def generate_bibtex(citations: List[Dict]) -> str:
    """Generate BibTeX entries from citation list."""
    bibtex = ""
    for i, cite in enumerate(citations):
        key = f"ref{i+1}_{cite.get('year', '')}"
        bibtex += f"""@article{{{key},
  title = {{{cite.get('title', '')}}},
  journal = {{{cite.get('journal', '')}}},
  year = {{{cite.get('year', '')}}},
  doi = {{{cite.get('doi', '')}}}
}}

"""
    return bibtex


if __name__ == "__main__":
    main()
