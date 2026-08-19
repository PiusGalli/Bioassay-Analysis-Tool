from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def draw_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    d = rdMolDraw2D.MolDraw2DSVG(300, 200)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return d.GetDrawingText()
