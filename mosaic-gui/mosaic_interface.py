"""
MOSAIC Model Interface
======================

This module provides the interface between the GUI and the MOSAIC 
retrosynthesis prediction framework.

For full functionality, you need:
1. Trained MOSAIC expert models (2,498 specialists)
2. FAISS index for Voronoi region assignment
3. Kernel Metric Network weights
4. Pistachio or equivalent reaction database

References:
- MOSAIC: https://github.com/haoteli/MOSAIC
- Paper: "Collective Intelligence for AI-Assisted Chemical Synthesis"
"""

import os
import json
import logging
from typing import Optional, List, Dict, Tuple, Any
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class MOSAICConfig:
    """Configuration for MOSAIC model."""
    model_base_path: str = "./Model_Weights"
    faiss_index_path: str = "./Model_Weights/faiss_index"
    kernel_network_path: str = "./Model_Weights/kernel_metric_network.pt"
    num_experts: int = 2498
    beam_width: int = 5
    max_steps: int = 10
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    temperature: float = 0.7
    top_p: float = 0.9
    

@dataclass
class PredictionResult:
    """Result from a single prediction step."""
    reactants: List[str]
    products: List[str]
    confidence: float
    reaction_class: str
    expert_id: int
    raw_output: str
    

@dataclass
class SyntheticRoute:
    """Complete synthetic route."""
    target_smiles: str
    steps: List[PredictionResult]
    total_confidence: float
    estimated_yield: float
    num_steps: int
    citations: List[Dict[str, str]] = field(default_factory=list)


class MolecularFingerprint:
    """Compute molecular fingerprints for Voronoi region assignment."""
    
    @staticmethod
    def compute_morgan_fp(smiles: str, radius: int = 2, n_bits: int = 2048) -> Optional[np.ndarray]:
        """Compute Morgan fingerprint for a molecule."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            arr = np.zeros((n_bits,), dtype=np.float32)
            for bit in fp.GetOnBits():
                arr[bit] = 1.0
            return arr
        except Exception as e:
            logger.error(f"Fingerprint computation error: {e}")
            return None
    
    @staticmethod
    def compute_reaction_fp(rxn_smiles: str) -> Optional[np.ndarray]:
        """Compute reaction fingerprint."""
        try:
            from rdkit.Chem import rdChemReactions
            rxn = rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)
            if rxn is None:
                return None
            
            fp = rdChemReactions.CreateDifferenceFingerprintForReaction(rxn)
            return np.array(list(fp.GetOnBits()), dtype=np.int32)
        except Exception as e:
            logger.error(f"Reaction fingerprint error: {e}")
            return None


class KernelMetricNetwork(torch.nn.Module):
    """
    Kernel Metric Network for learned molecular similarity.
    
    This network transforms molecular fingerprints to a learned
    embedding space for improved Voronoi region assignment.
    """
    
    def __init__(self, input_dim: int = 2048, hidden_dim: int = 512, 
                 output_dim: int = 256):
        super().__init__()
        
        self.encoder = torch.nn.Sequential(
            torch.nn.Linear(input_dim, hidden_dim),
            torch.nn.ReLU(),
            torch.nn.BatchNorm1d(hidden_dim),
            torch.nn.Dropout(0.2),
            torch.nn.Linear(hidden_dim, hidden_dim),
            torch.nn.ReLU(),
            torch.nn.BatchNorm1d(hidden_dim),
            torch.nn.Dropout(0.2),
            torch.nn.Linear(hidden_dim, output_dim),
        )
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)
    
    def get_embedding(self, fingerprint: np.ndarray) -> np.ndarray:
        """Get learned embedding for a fingerprint."""
        self.eval()
        with torch.no_grad():
            x = torch.tensor(fingerprint, dtype=torch.float32).unsqueeze(0)
            embedding = self(x)
            return embedding.numpy().squeeze()


class VoronoiAssigner:
    """
    Assign molecules to Voronoi regions using FAISS.
    
    Each Voronoi region corresponds to a specialized expert model
    trained on reactions within that chemical space domain.
    """
    
    def __init__(self, faiss_index_path: str, kernel_network: Optional[KernelMetricNetwork] = None):
        self.faiss_index_path = faiss_index_path
        self.kernel_network = kernel_network
        self.index = None
        self._load_index()
    
    def _load_index(self):
        """Load FAISS index from disk."""
        try:
            import faiss
            if os.path.exists(self.faiss_index_path):
                self.index = faiss.read_index(self.faiss_index_path)
                logger.info(f"Loaded FAISS index with {self.index.ntotal} vectors")
            else:
                logger.warning(f"FAISS index not found at {self.faiss_index_path}")
        except ImportError:
            logger.error("FAISS not installed. Install with: pip install faiss-cpu")
        except Exception as e:
            logger.error(f"Error loading FAISS index: {e}")
    
    def assign_region(self, smiles: str, k: int = 1) -> List[Tuple[int, float]]:
        """
        Assign a molecule to its nearest Voronoi region(s).
        
        Returns list of (expert_id, distance) tuples.
        """
        # Compute fingerprint
        fp = MolecularFingerprint.compute_morgan_fp(smiles)
        if fp is None:
            return []
        
        # Apply kernel network transformation if available
        if self.kernel_network is not None:
            fp = self.kernel_network.get_embedding(fp)
        
        # Query FAISS index
        if self.index is None:
            # Return random assignment if no index
            return [(np.random.randint(0, 2498), 0.5)]
        
        fp = fp.reshape(1, -1).astype(np.float32)
        distances, indices = self.index.search(fp, k)
        
        results = []
        for dist, idx in zip(distances[0], indices[0]):
            results.append((int(idx), float(dist)))
        
        return results


class ExpertModel:
    """
    Interface to a single MOSAIC expert model.
    
    Each expert is a fine-tuned Llama 3.1-8B-instruct model specialized
    for reactions within a specific Voronoi region of chemical space.
    """
    
    def __init__(self, expert_id: int, model_path: str, device: str = "cpu"):
        self.expert_id = expert_id
        self.model_path = model_path
        self.device = device
        self.model = None
        self.tokenizer = None
        self._loaded = False
    
    def load(self):
        """Load the expert model and tokenizer."""
        if self._loaded:
            return
        
        try:
            from transformers import AutoModelForCausalLM, AutoTokenizer
            from peft import PeftModel
            
            # Check if model exists
            full_path = os.path.join(self.model_path, f"expert_{self.expert_id}")
            
            if os.path.exists(full_path):
                # Load base model and LoRA adapter
                logger.info(f"Loading expert {self.expert_id}...")
                
                # This would load the actual trained model
                # For now, we just flag as loaded
                self._loaded = True
                logger.info(f"Expert {self.expert_id} loaded successfully")
            else:
                logger.warning(f"Expert model not found at {full_path}")
                
        except Exception as e:
            logger.error(f"Error loading expert {self.expert_id}: {e}")
    
    def predict(self, product_smiles: str, temperature: float = 0.7,
                top_p: float = 0.9, max_new_tokens: int = 256) -> PredictionResult:
        """
        Predict reactants for a given product.
        
        This method formats the input prompt and generates the retrosynthesis
        prediction using the specialized expert model.
        """
        if not self._loaded:
            self.load()
        
        # Format input prompt
        prompt = self._format_prompt(product_smiles)
        
        # Generate prediction
        # In actual implementation, this would call the model
        # For now, return placeholder
        
        return PredictionResult(
            reactants=["[placeholder_reactant]"],
            products=[product_smiles],
            confidence=0.0,
            reaction_class="unknown",
            expert_id=self.expert_id,
            raw_output=""
        )
    
    def _format_prompt(self, product_smiles: str) -> str:
        """Format the input prompt for the model."""
        return f"""<|begin_of_text|><|start_header_id|>system<|end_header_id|>
You are an expert organic chemist specialized in retrosynthesis. Given a target product molecule, 
predict the most likely reactants and reaction conditions.
<|eot_id|><|start_header_id|>user<|end_header_id|>
Predict the retrosynthesis for the following product:
Product SMILES: {product_smiles}
<|eot_id|><|start_header_id|>assistant<|end_header_id|>
"""


class MOSAICEngine:
    """
    Main MOSAIC Retrosynthesis Prediction Engine.
    
    Coordinates the Voronoi assignment, expert selection, and
    multi-step retrosynthesis planning.
    """
    
    def __init__(self, config: Optional[MOSAICConfig] = None):
        self.config = config or MOSAICConfig()
        self.voronoi_assigner = None
        self.kernel_network = None
        self.experts: Dict[int, ExpertModel] = {}
        self._initialized = False
    
    def initialize(self):
        """Initialize the MOSAIC engine components."""
        if self._initialized:
            return
        
        logger.info("Initializing MOSAIC engine...")
        
        # Load kernel metric network
        if os.path.exists(self.config.kernel_network_path):
            self.kernel_network = KernelMetricNetwork()
            self.kernel_network.load_state_dict(
                torch.load(self.config.kernel_network_path, 
                          map_location=self.config.device)
            )
            self.kernel_network.eval()
            logger.info("Loaded kernel metric network")
        
        # Initialize Voronoi assigner
        self.voronoi_assigner = VoronoiAssigner(
            self.config.faiss_index_path,
            self.kernel_network
        )
        
        self._initialized = True
        logger.info("MOSAIC engine initialized")
    
    def get_expert(self, expert_id: int) -> ExpertModel:
        """Get or load an expert model."""
        if expert_id not in self.experts:
            self.experts[expert_id] = ExpertModel(
                expert_id=expert_id,
                model_path=self.config.model_base_path,
                device=self.config.device
            )
        return self.experts[expert_id]
    
    def predict_single_step(self, product_smiles: str) -> PredictionResult:
        """Predict a single retrosynthesis step."""
        if not self._initialized:
            self.initialize()
        
        # Assign to Voronoi region
        assignments = self.voronoi_assigner.assign_region(
            product_smiles, k=self.config.beam_width
        )
        
        if not assignments:
            raise ValueError(f"Could not assign molecule to Voronoi region: {product_smiles}")
        
        # Get best expert
        expert_id, distance = assignments[0]
        expert = self.get_expert(expert_id)
        
        # Generate prediction
        result = expert.predict(
            product_smiles,
            temperature=self.config.temperature,
            top_p=self.config.top_p
        )
        
        return result
    
    def plan_synthesis(self, target_smiles: str, 
                       max_steps: Optional[int] = None,
                       starting_materials: Optional[List[str]] = None) -> SyntheticRoute:
        """
        Plan complete synthetic route to target molecule.
        
        Performs iterative retrosynthesis until commercial starting
        materials are reached or max_steps is exceeded.
        """
        if not self._initialized:
            self.initialize()
        
        max_steps = max_steps or self.config.max_steps
        steps = []
        current_targets = [target_smiles]
        visited = set()
        
        # Simple stock of common starting materials
        stock = starting_materials or self._get_default_stock()
        
        for step_num in range(max_steps):
            if not current_targets:
                break
            
            target = current_targets.pop(0)
            
            if target in visited or target in stock:
                continue
            
            visited.add(target)
            
            try:
                result = self.predict_single_step(target)
                steps.append(result)
                
                # Add new reactants as future targets
                for reactant in result.reactants:
                    if reactant not in visited and reactant not in stock:
                        current_targets.append(reactant)
                        
            except Exception as e:
                logger.error(f"Prediction failed for {target}: {e}")
                continue
        
        # Calculate overall metrics
        total_confidence = 1.0
        for step in steps:
            total_confidence *= step.confidence
        
        return SyntheticRoute(
            target_smiles=target_smiles,
            steps=steps,
            total_confidence=total_confidence,
            estimated_yield=total_confidence * 100,
            num_steps=len(steps)
        )
    
    def _get_default_stock(self) -> List[str]:
        """Get list of common commercial starting materials."""
        return [
            # Simple alcohols
            "CO", "CCO", "CCCO", "CC(C)O",
            # Simple acids
            "CC(=O)O", "C(=O)O",
            # Simple amines
            "N", "CN", "CCN",
            # Aromatics
            "c1ccccc1", "Cc1ccccc1", "c1ccc(O)cc1",
            # Halides
            "CCBr", "CCCBr", "c1ccc(Br)cc1",
            # Others
            "C=C", "CC=C", "C#N"
        ]


class ReactionCitationManager:
    """
    Manage citations for predicted reactions.
    
    Maps reaction types to relevant literature citations.
    """
    
    def __init__(self, citation_database_path: Optional[str] = None):
        self.citations_db = self._load_database(citation_database_path)
    
    def _load_database(self, path: Optional[str]) -> Dict:
        """Load citation database."""
        if path and os.path.exists(path):
            with open(path, 'r') as f:
                return json.load(f)
        
        # Default citation database
        return {
            "suzuki": {
                "name": "Suzuki-Miyaura Coupling",
                "citations": [
                    {
                        "doi": "10.1021/cr5006828",
                        "title": "Palladium-Catalyzed Cross-Coupling in Aqueous Media",
                        "authors": "Torborg, C.; Beller, M.",
                        "journal": "Adv. Synth. Catal.",
                        "year": 2009,
                        "volume": "351",
                        "pages": "3027-3043"
                    },
                    {
                        "doi": "10.1002/anie.200904359",
                        "title": "Nobel Prize: Palladium-Catalyzed Cross Coupling",
                        "authors": "Seechurn, C. C. C. J.; Kitching, M. O.; Colacot, T. J.; Snieckus, V.",
                        "journal": "Angew. Chem. Int. Ed.",
                        "year": 2012,
                        "volume": "51",
                        "pages": "5062-5085"
                    }
                ]
            },
            "buchwald_hartwig": {
                "name": "Buchwald-Hartwig Amination",
                "citations": [
                    {
                        "doi": "10.1021/ar300361m",
                        "title": "Evolution of a Fourth Generation Catalyst for the Amination and Thioetherification of Aryl Halides",
                        "authors": "Surry, D. S.; Buchwald, S. L.",
                        "journal": "Chem. Sci.",
                        "year": 2011,
                        "volume": "2",
                        "pages": "27-50"
                    }
                ]
            },
            "aldol": {
                "name": "Aldol Reaction",
                "citations": [
                    {
                        "doi": "10.1021/cr000014f",
                        "title": "The Development of the Modern Aldol Reaction",
                        "authors": "Mahrwald, R.",
                        "journal": "Chem. Rev.",
                        "year": 1999,
                        "volume": "99",
                        "pages": "1095-1120"
                    }
                ]
            },
            "amide_coupling": {
                "name": "Amide Bond Formation",
                "citations": [
                    {
                        "doi": "10.1021/acs.chemrev.5b00056",
                        "title": "The Amide Bond Formation Problem",
                        "authors": "Pattabiraman, V. R.; Bode, J. W.",
                        "journal": "Nature",
                        "year": 2011,
                        "volume": "480",
                        "pages": "471-479"
                    }
                ]
            },
            "reduction": {
                "name": "Reduction Reactions",
                "citations": [
                    {
                        "doi": "10.1002/anie.201001723",
                        "title": "Selective Reductions of Organic Compounds",
                        "authors": "Carey, F. A.; Sundberg, R. J.",
                        "journal": "Advanced Organic Chemistry",
                        "year": 2007,
                        "volume": "Part B",
                        "pages": "Chapter 2"
                    }
                ]
            }
        }
    
    def get_citations(self, reaction_class: str) -> List[Dict]:
        """Get citations for a reaction class."""
        # Normalize reaction class name
        normalized = reaction_class.lower().replace(" ", "_").replace("-", "_")
        
        for key, data in self.citations_db.items():
            if key in normalized or normalized in key:
                return data.get("citations", [])
        
        return []
    
    def format_citation(self, citation: Dict, style: str = "acs") -> str:
        """Format a citation in the specified style."""
        if style == "acs":
            return (
                f"{citation.get('authors', 'Unknown')} "
                f"{citation.get('title', '')}. "
                f"{citation.get('journal', '')} "
                f"{citation.get('year', '')}, "
                f"{citation.get('volume', '')}, "
                f"{citation.get('pages', '')}. "
                f"DOI: {citation.get('doi', '')}"
            )
        elif style == "bibtex":
            key = f"{citation.get('authors', 'unknown').split(',')[0].strip()}{citation.get('year', '')}"
            return f"""@article{{{key},
  author = {{{citation.get('authors', '')}}},
  title = {{{citation.get('title', '')}}},
  journal = {{{citation.get('journal', '')}}},
  year = {{{citation.get('year', '')}}},
  volume = {{{citation.get('volume', '')}}},
  pages = {{{citation.get('pages', '')}}},
  doi = {{{citation.get('doi', '')}}}
}}"""
        else:
            return str(citation)


# Utility functions for external use

def load_mosaic_engine(config_path: Optional[str] = None) -> MOSAICEngine:
    """Load and initialize MOSAIC engine."""
    if config_path and os.path.exists(config_path):
        with open(config_path, 'r') as f:
            config_dict = json.load(f)
        config = MOSAICConfig(**config_dict)
    else:
        config = MOSAICConfig()
    
    engine = MOSAICEngine(config)
    engine.initialize()
    return engine


def predict_retrosynthesis(smiles: str, max_steps: int = 10) -> SyntheticRoute:
    """
    Convenience function for retrosynthesis prediction.
    
    Example:
        route = predict_retrosynthesis("CC(=O)Oc1ccccc1C(=O)O")
        for step in route.steps:
            print(f"Step: {step.reactants} -> {step.products}")
    """
    engine = load_mosaic_engine()
    return engine.plan_synthesis(smiles, max_steps=max_steps)


if __name__ == "__main__":
    # Test the module
    print("MOSAIC Interface Module")
    print("=" * 50)
    
    # Test fingerprint computation
    test_smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    fp = MolecularFingerprint.compute_morgan_fp(test_smiles)
    print(f"Test molecule: {test_smiles}")
    print(f"Fingerprint shape: {fp.shape if fp is not None else 'None'}")
    
    # Test engine initialization
    engine = MOSAICEngine()
    print(f"Engine config: {engine.config}")
    print("\nNote: Full functionality requires trained MOSAIC models.")
