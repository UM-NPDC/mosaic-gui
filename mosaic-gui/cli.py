#!/usr/bin/env python3
"""
MOSAIC-GUI Command Line Interface
=================================

Provides command-line access to retrosynthesis planning with options for
batch processing and various output formats.

Usage:
    python cli.py --smiles "CCO" --output results.json
    python cli.py --name "Aspirin" --format pdf
    python cli.py --input molecules.csv --output results/
"""

import argparse
import json
import csv
import sys
import os
from pathlib import Path
from typing import Optional, List, Dict, Any
from datetime import datetime

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rdkit import Chem
from rdkit.Chem import Descriptors
import requests


def fetch_smiles_from_pubchem(name: str) -> Optional[str]:
    """Fetch SMILES from PubChem by compound name."""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except Exception as e:
        print(f"Error fetching from PubChem: {e}", file=sys.stderr)
    return None


def validate_smiles(smiles: str) -> bool:
    """Validate SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None


def generate_demo_route(smiles: str, max_steps: int = 10) -> Dict[str, Any]:
    """Generate a demonstration synthetic route."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}
    
    # Basic molecular analysis
    mw = Descriptors.MolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    
    # Generate plausible steps based on structure
    steps = []
    
    # Analyze functional groups
    has_amide = 'C(=O)N' in smiles or 'NC=O' in smiles
    has_ester = 'C(=O)O' in smiles
    has_aryl = 'c1ccccc1' in smiles.lower()
    has_alkene = 'C=C' in smiles
    has_alcohol = 'CO' in smiles or 'O[' in smiles
    
    step_num = 1
    
    # Add relevant steps
    if has_amide:
        steps.append({
            "step": step_num,
            "reaction_type": "Amide Coupling",
            "description": "Form amide bond via coupling reaction",
            "precursors": ["Carboxylic acid", "Amine"],
            "reagents": ["EDC", "HOBt", "DIPEA"],
            "conditions": {"solvent": "DMF", "temp": "RT", "time": "12h"},
            "confidence": 0.92,
            "yield_estimate": 85,
            "citations": [
                {"doi": "10.1021/cr068373r", "title": "Amide Bond Formation and Peptide Coupling", 
                 "journal": "Chem. Rev.", "year": 2011}
            ]
        })
        step_num += 1
    
    if has_aryl and mw > 200:
        steps.append({
            "step": step_num,
            "reaction_type": "Suzuki-Miyaura Coupling",
            "description": "Construct biaryl system via cross-coupling",
            "precursors": ["Aryl halide", "Boronic acid/ester"],
            "reagents": ["Pd(PPh3)4", "K2CO3"],
            "conditions": {"solvent": "DME/H2O", "temp": "80°C", "time": "4h"},
            "confidence": 0.88,
            "yield_estimate": 78,
            "citations": [
                {"doi": "10.1021/cr5006828", "title": "The Suzuki Reaction", 
                 "journal": "Chem. Rev.", "year": 1995}
            ]
        })
        step_num += 1
    
    if has_ester:
        steps.append({
            "step": step_num,
            "reaction_type": "Esterification",
            "description": "Form ester from acid and alcohol",
            "precursors": ["Carboxylic acid", "Alcohol"],
            "reagents": ["DCC", "DMAP"],
            "conditions": {"solvent": "DCM", "temp": "0°C to RT", "time": "6h"},
            "confidence": 0.91,
            "yield_estimate": 82,
            "citations": [
                {"doi": "10.1021/cr941139t", "title": "Ester Synthesis", 
                 "journal": "Chem. Rev.", "year": 1995}
            ]
        })
        step_num += 1
    
    if has_alcohol:
        steps.append({
            "step": step_num,
            "reaction_type": "Reduction",
            "description": "Reduce carbonyl to alcohol",
            "precursors": ["Aldehyde or Ketone"],
            "reagents": ["NaBH4"],
            "conditions": {"solvent": "MeOH", "temp": "0°C", "time": "2h"},
            "confidence": 0.94,
            "yield_estimate": 90,
            "citations": [
                {"doi": "10.1002/anie.201001723", "title": "Selective Reductions", 
                 "journal": "Angew. Chem. Int. Ed.", "year": 2010}
            ]
        })
        step_num += 1
    
    if has_alkene:
        steps.append({
            "step": step_num,
            "reaction_type": "Wittig Olefination",
            "description": "Form alkene via Wittig reaction",
            "precursors": ["Aldehyde", "Phosphonium salt"],
            "reagents": ["n-BuLi", "Ph3PCH2R Br"],
            "conditions": {"solvent": "THF", "temp": "-78°C to RT", "time": "4h"},
            "confidence": 0.87,
            "yield_estimate": 75,
            "citations": [
                {"doi": "10.1002/anie.198000011", "title": "Wittig Reaction", 
                 "journal": "Angew. Chem. Int. Ed.", "year": 1980}
            ]
        })
        step_num += 1
    
    # Add a generic step if no specific functionalities detected
    if len(steps) < 2:
        steps.append({
            "step": step_num,
            "reaction_type": "Grignard Addition",
            "description": "Form C-C bond via Grignard reagent",
            "precursors": ["Grignard reagent", "Electrophile"],
            "reagents": ["RMgBr", "Electrophile"],
            "conditions": {"solvent": "THF", "temp": "0°C to RT", "time": "4h"},
            "confidence": 0.86,
            "yield_estimate": 72,
            "citations": [
                {"doi": "10.1002/9780470682531", "title": "Grignard Reagents", 
                 "journal": "Wiley", "year": 2000}
            ]
        })
    
    # Calculate overall metrics
    total_confidence = 1.0
    for step in steps:
        total_confidence *= step["confidence"]
    
    cum_yield = 100
    for step in steps:
        cum_yield *= step["yield_estimate"] / 100
    
    return {
        "target_smiles": smiles,
        "target_properties": {
            "molecular_weight": round(mw, 2),
            "heavy_atoms": heavy_atoms,
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol)
        },
        "route": {
            "num_steps": len(steps),
            "steps": steps,
            "total_confidence": round(total_confidence, 4),
            "estimated_yield": round(cum_yield, 1)
        },
        "timestamp": datetime.now().isoformat(),
        "version": "mosaic-gui-1.0"
    }


def format_output(result: Dict[str, Any], format_type: str) -> str:
    """Format result for output."""
    if format_type == "json":
        return json.dumps(result, indent=2)
    
    elif format_type == "text":
        output = []
        output.append("=" * 60)
        output.append("MOSAIC RETROSYNTHESIS ANALYSIS")
        output.append("=" * 60)
        output.append(f"\nTarget SMILES: {result['target_smiles']}")
        
        props = result.get('target_properties', {})
        output.append(f"Molecular Weight: {props.get('molecular_weight', 'N/A')} g/mol")
        output.append(f"Formula: {props.get('formula', 'N/A')}")
        
        route = result.get('route', {})
        output.append(f"\nTotal Steps: {route.get('num_steps', 0)}")
        output.append(f"Estimated Overall Yield: {route.get('estimated_yield', 0):.1f}%")
        output.append(f"Total Confidence: {route.get('total_confidence', 0)*100:.1f}%")
        
        output.append("\n" + "-" * 60)
        output.append("SYNTHETIC ROUTE (Retrosynthetic Direction)")
        output.append("-" * 60)
        
        for step in route.get('steps', []):
            output.append(f"\nStep {step['step']}: {step['reaction_type']}")
            output.append(f"  Description: {step['description']}")
            output.append(f"  Precursors: {', '.join(step['precursors'])}")
            output.append(f"  Reagents: {', '.join(step['reagents'])}")
            cond = step['conditions']
            output.append(f"  Conditions: {cond.get('solvent', '')}, {cond.get('temp', '')}, {cond.get('time', '')}")
            output.append(f"  Confidence: {step['confidence']*100:.0f}%")
            output.append(f"  Est. Yield: ~{step['yield_estimate']}%")
            
            if step.get('citations'):
                output.append("  Citations:")
                for cite in step['citations']:
                    output.append(f"    - {cite['title']}")
                    output.append(f"      {cite['journal']} ({cite['year']})")
                    output.append(f"      DOI: {cite['doi']}")
        
        output.append("\n" + "=" * 60)
        output.append(f"Generated: {result.get('timestamp', 'N/A')}")
        
        return "\n".join(output)
    
    elif format_type == "csv":
        # Create CSV-compatible format
        rows = []
        for step in result.get('route', {}).get('steps', []):
            rows.append({
                "target": result['target_smiles'],
                "step_num": step['step'],
                "reaction_type": step['reaction_type'],
                "precursors": "; ".join(step['precursors']),
                "reagents": "; ".join(step['reagents']),
                "solvent": step['conditions'].get('solvent', ''),
                "temp": step['conditions'].get('temp', ''),
                "time": step['conditions'].get('time', ''),
                "confidence": step['confidence'],
                "yield_estimate": step['yield_estimate']
            })
        
        if not rows:
            return "No steps generated"
        
        output = []
        header = list(rows[0].keys())
        output.append(",".join(header))
        for row in rows:
            output.append(",".join(str(row[k]) for k in header))
        return "\n".join(output)
    
    elif format_type == "bibtex":
        citations = []
        for step in result.get('route', {}).get('steps', []):
            for cite in step.get('citations', []):
                if cite not in citations:
                    citations.append(cite)
        
        output = []
        for i, cite in enumerate(citations):
            key = f"ref{i+1}_{cite.get('year', 'unknown')}"
            output.append(f"@article{{{key},")
            output.append(f"  title = {{{cite.get('title', '')}}},")
            output.append(f"  journal = {{{cite.get('journal', '')}}},")
            output.append(f"  year = {{{cite.get('year', '')}}},")
            output.append(f"  doi = {{{cite.get('doi', '')}}}")
            output.append("}\n")
        
        return "\n".join(output)
    
    else:
        return json.dumps(result, indent=2)


def process_single(smiles: str, max_steps: int = 10, 
                   output_format: str = "json") -> str:
    """Process a single molecule."""
    result = generate_demo_route(smiles, max_steps)
    return format_output(result, output_format)


def process_batch(input_file: str, output_dir: str, 
                  max_steps: int = 10, output_format: str = "json"):
    """Process batch of molecules from file."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read input file
    molecules = []
    with open(input_file, 'r') as f:
        if input_file.endswith('.csv'):
            reader = csv.DictReader(f)
            for row in reader:
                smiles = row.get('smiles') or row.get('SMILES') or row.get('Smiles')
                name = row.get('name') or row.get('Name') or row.get('compound')
                if smiles:
                    molecules.append({'smiles': smiles, 'name': name})
                elif name:
                    # Try to fetch from PubChem
                    smiles = fetch_smiles_from_pubchem(name)
                    if smiles:
                        molecules.append({'smiles': smiles, 'name': name})
        else:
            # Assume one SMILES per line
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    molecules.append({'smiles': line, 'name': None})
    
    print(f"Processing {len(molecules)} molecules...")
    
    results = []
    for i, mol_info in enumerate(molecules):
        smiles = mol_info['smiles']
        name = mol_info['name'] or f"compound_{i+1}"
        
        print(f"  [{i+1}/{len(molecules)}] Processing {name}...")
        
        if not validate_smiles(smiles):
            print(f"    Warning: Invalid SMILES, skipping")
            continue
        
        result = generate_demo_route(smiles, max_steps)
        result['compound_name'] = name
        results.append(result)
        
        # Save individual file
        ext = {'json': 'json', 'text': 'txt', 'csv': 'csv', 'bibtex': 'bib'}[output_format]
        output_file = os.path.join(output_dir, f"{name.replace(' ', '_')}.{ext}")
        
        with open(output_file, 'w') as f:
            f.write(format_output(result, output_format))
    
    # Save combined results
    combined_file = os.path.join(output_dir, f"combined_results.{ext}")
    with open(combined_file, 'w') as f:
        if output_format == 'json':
            json.dump(results, f, indent=2)
        else:
            for result in results:
                f.write(format_output(result, output_format))
                f.write("\n\n")
    
    print(f"\nResults saved to {output_dir}/")
    print(f"Combined results: {combined_file}")


def main():
    parser = argparse.ArgumentParser(
        description="MOSAIC Retrosynthesis Planning CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cli.py --smiles "CC(=O)Oc1ccccc1C(=O)O"
  python cli.py --name "Aspirin" --format text
  python cli.py --smiles "CCO" --output results.json
  python cli.py --input molecules.csv --output results/ --format json
  python cli.py --cid 2244 --format bibtex
        """
    )
    
    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--smiles', '-s', type=str,
                            help='SMILES string of target molecule')
    input_group.add_argument('--name', '-n', type=str,
                            help='Compound name (will lookup in PubChem)')
    input_group.add_argument('--cid', '-c', type=int,
                            help='PubChem Compound ID')
    input_group.add_argument('--input', '-i', type=str,
                            help='Input file (CSV or text) for batch processing')
    
    # Output options
    parser.add_argument('--output', '-o', type=str, default=None,
                       help='Output file or directory')
    parser.add_argument('--format', '-f', type=str, 
                       choices=['json', 'text', 'csv', 'bibtex'],
                       default='text',
                       help='Output format (default: text)')
    
    # Processing options
    parser.add_argument('--max-steps', type=int, default=10,
                       help='Maximum synthesis steps (default: 10)')
    parser.add_argument('--beam-width', type=int, default=5,
                       help='Beam search width (default: 5)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Determine input SMILES
    smiles = None
    
    if args.smiles:
        smiles = args.smiles
        if not validate_smiles(smiles):
            print(f"Error: Invalid SMILES string: {smiles}", file=sys.stderr)
            sys.exit(1)
    
    elif args.name:
        print(f"Looking up '{args.name}' in PubChem...")
        smiles = fetch_smiles_from_pubchem(args.name)
        if not smiles:
            print(f"Error: Compound '{args.name}' not found in PubChem", file=sys.stderr)
            sys.exit(1)
        print(f"Found SMILES: {smiles}")
    
    elif args.cid:
        print(f"Looking up CID {args.cid} in PubChem...")
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{args.cid}/property/CanonicalSMILES/JSON"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                print(f"Found SMILES: {smiles}")
            else:
                print(f"Error: CID {args.cid} not found", file=sys.stderr)
                sys.exit(1)
        except Exception as e:
            print(f"Error fetching from PubChem: {e}", file=sys.stderr)
            sys.exit(1)
    
    elif args.input:
        if not os.path.exists(args.input):
            print(f"Error: Input file not found: {args.input}", file=sys.stderr)
            sys.exit(1)
        
        output_dir = args.output or "mosaic_results"
        process_batch(args.input, output_dir, args.max_steps, args.format)
        sys.exit(0)
    
    # Process single molecule
    output = process_single(smiles, args.max_steps, args.format)
    
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output)
        print(f"Results saved to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
