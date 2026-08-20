from rdkit import Chem, DataStructs
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdFingerprintGenerator
import numpy as np
import pandas as pd

def draw_smiles(smiles):
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
      return None
  d = rdMolDraw2D.MolDraw2DSVG(300, 200)
  d.DrawMolecule(mol)
  d.FinishDrawing()
  return d.GetDrawingText()

def smiles_to_mols(smiles_list, labels):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    
    valid_mols = [m for m in mols if m is not None]
    valid_indices = [i for i, m in enumerate(mols) if m is not None]
    
    valid_labels = [labels[i] for i in valid_indices]
    
    none_indices = [i for i, m in enumerate(mols) if m is None]
    
    return valid_mols, none_indices, valid_labels

def morgan_fingerprint(mols, radius=2, nBits=1024):
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=nBits,
        includeChirality=True,
    )
    mol_fingerprints = [generator.GetFingerprint(mol) for mol in mols]
    return mol_fingerprints

def tanimoto_similarity(mol_fingerprints, labels):
  n = len(mol_fingerprints)
  tanimoto_similarities = np.ones((n,n))
  
  for i in range(len(mol_fingerprints)):
    for j in range(i+1, len(mol_fingerprints)):
      mol1 = mol_fingerprints[i]
      mol2 = mol_fingerprints[j]
      
      tanimoto = DataStructs.TanimotoSimilarity(mol1, mol2)
      tanimoto_similarities[i,j] = tanimoto
      tanimoto_similarities[j,i] = tanimoto
  
  return pd.DataFrame(tanimoto_similarities, index=labels, columns=labels)

def run_similarity_pipeline(smiles_list, labels, radius=2, nBits=1024):
    valid_mols, failed_indices, valid_labels = smiles_to_mols(smiles_list, labels)
    fingerprints = morgan_fingerprint(valid_mols, radius=radius, nBits=nBits)
    sim_matrix = tanimoto_similarity(fingerprints, valid_labels)
    
    return {
        "similarity_matrix": sim_matrix,
        "failed_indices": failed_indices,
        "n_valid": len(valid_mols)
    }

