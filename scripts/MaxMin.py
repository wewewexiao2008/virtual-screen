from rdkit import Chem
from rdkit.Chem import AllChem as ch
from rdkit.Chem import Draw as d
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

mols = ch.SDMolSupplier('/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/Alinda/Alinda_StockSC_2021_December_202420_sdf.sdf')
mol1 = [x for x in mols if x is not None]
print('Number:' + str(len(mol1)))
fps = [ch.GetMorganFingerprint(m, 3) for m in mol1]

def distij(i, j, fps=fps):
    return 1 - DataStructs.DiceSimilarity(fps[i], fps[j])
picker = MaxMinPicker()
pickIndices = picker.LazyPick(distij, len(fps), 10, seed=23)
picks = [mol1[x] for x in pickIndices]
len(picks)
