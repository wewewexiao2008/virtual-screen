#!/usr/bin/env python
#
# 
from rdkit import Chem
from rdkit.Chem import AllChem as ch
from rdkit.Chem import Draw as d
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs


#PATH and NAME of SDF
PATH1 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/Alinda/'
NAME1 = 'Alinda'
PATH2 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/Analytion'
NAME2 = 'Analytion'
PATH3 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/Aronis'
NAME3 = 'Aronis'
PATH4 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/Asinex'
NAME4 = 'Asinex'
PATH5 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/BIONET'
NAME5 = 'BIONET'
PATH6 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/Chembridge'
NAME6 = 'Chembridge'
PATH7 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/ChemicalBlock'
NAME7 = 'ChemicalBlock'
PATH8 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/D009'
NAME8 = 'D009'
PATH9 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/D063'
NAME9 = 'D063'
PATH10 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/IBS'
NAME10 = 'IBS'
PATH11 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/known'
NAME11 = 'known'
PATH12 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/Lifechemicals'
NAME12 = 'Lifechemicals'
PATH13 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/Maybridge'
NAME13 = 'Maybridge'
PATH14 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/Princeton'
NAME14 = 'Princeton'
PATH15 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/Specs'
NAME15 = 'Specs'
PATH16 = '/work/home/Uranus/data/Virtualflow/VFLP/ligand-library/fingerprint/SDF/T001'
NAME16 = 'T001'

#Read SDF
globals()[str(NAME1)] = ch.SDMolSupplier('%s/Alinda_StockSC_2021_December_202420_sdf.sdf' %PATH1)
mol1 = [x for x in globals()[str(NAME1)] if x is not None]
print('Number of %s: ' %NAME1 + str(len(mol1)))

globals()[str(NAME2) + "_1"] = ch.SDMolSupplier('%s/MACROx_4587.sdf' %PATH2)
globals()[str(NAME2) + "_2"] = ch.SDMolSupplier('%s/MEGx_414.sdf' %PATH2)
globals()[str(NAME2) + "_3"] = ch.SDMolSupplier('%s/MEGx_6348.sdf' %PATH2)
globals()[str(NAME2) + "_4"] = ch.SDMolSupplier('%s/NATx_33222.sdf' %PATH2)
globals()[str(NAME2) + "_5"] = ch.SDMolSupplier('%s/NATx_965.sdf' %PATH2)
mol2 = [x for x in globals()[str(NAME2) + "_1"] if x is not None] + [x for x in globals()[str(NAME2) + "_2"] if x is not None] + [x for x in globals()[str(NAME2) + "_3"] if x is not None] + [x for x in globals()[str(NAME2) + "_4"] if x is not None] + [x for x in globals()[str(NAME2) + "_5"] if x is not None]
print('Number of %s: ' %NAME2 + str(len(mol2)))

globals()[str(NAME3)] = ch.SDMolSupplier('%s/ScreeningStock.sdf' %PATH3)
mol3 = [x for x in globals()[str(NAME3)] if x is not None]
print('Number of %s: ' %NAME3 + str(len(mol3)))s

globals()[str(NAME4)] = ch.SDMolSupplier('%s/asinex_03_Feb_2022.sdf' %PATH4)
mol4 = [x for x in globals()[str(NAME4)] if x is not None]
print('Number of %s: ' %NAME4 + str(len(mol4)))

globals()[str(NAME5)] = ch.SDMolSupplier('%s/BIONET-February2022.sdf' %PATH5)
mol5 = [x for x in globals()[str(NAME5)] if x is not None]
print('Number of %s: ' %NAME5 + str(len(mol5)))

globals()[str(NAME6)] = ch.SDMolSupplier('%s/SC_Collection_Dec-2021.sdf' %PATH6)
mol6 = [x for x in globals()[str(NAME6)] if x is not None]
print('Number of %s: ' %NAME6 + str(len(mol6)))

globals()[str(NAME7)] = ch.SDMolSupplier('%s/screening_march2015.sdf' %PATH7)
mol7 = [x for x in globals()[str(NAME7)] if x is not None]
print('Number of %s: ' %NAME7 + str(len(mol7)))

globals()[str(NAME8) + "_1"] = ch.SDMolSupplier('%s/D009-1.sdf' %PATH8)
globals()[str(NAME8) + "_2"] = ch.SDMolSupplier('%s/D009-2.sdf' %PATH8)
globals()[str(NAME8) + "_3"] = ch.SDMolSupplier('%s/D009-3.sdf' %PATH8)
mol8 = [x for x in globals()[str(NAME8) + "_1"] if x is not None] + [x for x in globals()[str(NAME8) + "_2"] if x is not None] + [x for x in globals()[str(NAME8) + "_3"] if x is not None]
print('Number of %s: ' %NAME8 + str(len(mol8)))

globals()[str(NAME9)] = ch.SDMolSupplier('%s/D063.sdf' %PATH9)
mol9 = [x for x in globals()[str(NAME9)] if x is not None]
print('Number of %s: ' %NAME9 + str(len(mol9)))

globals()[str(NAME10) + "_1"] = ch.SDMolSupplier('%s/ibs2022feb_sc1.sdf' %PATH10)
globals()[str(NAME10) + "_2"] = ch.SDMolSupplier('%s/ibs2022feb_sc2.sdf' %PATH10)
globals()[str(NAME10) + "_3"] = ch.SDMolSupplier('%s/ibs2022feb_sc3.sdf' %PATH10)
mol10 = [x for x in globals()[str(NAME10) + "_1"] if x is not None] + [x for x in globals()[str(NAME10) + "_2"] if x is not None] + [x for x in globals()[str(NAME10) + "_3"] if x is not None]
print('Number of %s: ' %NAME10 + str(len(mol10)))

globals()[str(NAME11) + "_1"] = ch.SDMolSupplier('%s/D7800-Bioactive.sdf' %PATH11)
globals()[str(NAME11) + "_2"] = ch.SDMolSupplier('%s/L1000-TargetmolApproved.sdf' %PATH11)
globals()[str(NAME11) + "_3"] = ch.SDMolSupplier('%s/L4000-TargetmolBioactive.sdf' %PATH11)
mol11 = [x for x in globals()[str(NAME11) + "_1"] if x is not None] + [x for x in globals()[str(NAME11) + "_2"] if x is not None] + [x for x in globals()[str(NAME11) + "_3"] if x is not None]
print('Number of %s: ' %NAME11 + str(len(mol11)))

globals()[str(NAME12)] = ch.SDMolSupplier('%s/LC_Stock_HTS_Compounds.sdf' %PATH12)
mol12 = [x for x in globals()[str(NAME12)] if x is not None]
print('Number of %s: ' %NAME12 + str(len(mol12)))

globals()[str(NAME13)] = ch.SDMolSupplier('%s/Maybridge.sdf' %PATH13)
mol13 = [x for x in globals()[str(NAME13)] if x is not None]
print('Number of %s: ' %NAME13 + str(len(mol13)))

globals()[str(NAME14) + "_1"] = ch.SDMolSupplier('%s/Princeton_Screeneng_032020-001.sdf' %PATH14)
globals()[str(NAME14) + "_2"] = ch.SDMolSupplier('%s/Princeton_Screeneng_032020-002.sdf' %PATH14)
globals()[str(NAME14) + "_3"] = ch.SDMolSupplier('%s/Princeton_Screeneng_032020-003.sdf' %PATH14)
mol14 = [x for x in globals()[str(NAME14) + "_1"] if x is not None] + [x for x in globals()[str(NAME14) + "_2"] if x is not None] + [x for x in globals()[str(NAME14) + "_3"] if x is not None]
print('Number of %s: ' %NAME14 + str(len(mol14)))

globals()[str(NAME15)] = ch.SDMolSupplier('%s/Specs_SC_10mg_Mar2022.sdf' %PATH15)
mol15 = [x for x in globals()[str(NAME15)] if x is not None]
print('Number of %s: ' %NAME15 + str(len(mol15)))

globals()[str(NAME16) + "_1"] = ch.SDMolSupplier('%s/T001-1.sdf' %PATH16)
globals()[str(NAME16) + "_2"] = ch.SDMolSupplier('%s/T001-2.sdf' %PATH16)
mol16 = [x for x in globals()[str(NAME16) + "_1"] if x is not None] + [x for x in globals()[str(NAME16) + "_2"] if x is not None]
print('Number of %s: ' %NAME16 + str(len(mol16)))

#merge and calculate fingerprints of molecules
mols = mol1 + mol2 + mol3 + mol4 + mol5 + mol6 + mol7 + mol8 + mol9 + mol10 + mol11 + mol12 + mol13 + mol14 + mol15 + mol16
print()
print('SUM: ' + str(len(mols)))
print()
print('------------------------------------------------------------')
mols_fps_Morgen = [(m, ch.GetMorganFingerprintAsBitVect(m, 2)) for m in mols]

#type the smile of target molecule and calculate fingerprint
target = ch.MolFromSmiles('C1CN(C)CCN1CCCNC(=O)c2cc(=O)n(CC)c(c23)n(C)c(=O)n(c3=O)CC')
target_Morgen = ch.GetMorganFingerprintAsBitVect(target, 2)
target_MACCS = MACCSkeys.GenMACCSKeys(target)

#calculate similarity between target and ligand-library, extract top N molecules
print('Calculating Similarity...')
mols_targetsim_Morgen = [(m, DataStructs.TanimotoSimilarity(fp, target_Morgen)) for m, fp in mols_fps_Morgen]
sorted_mols_targetsim_Morgen = sorted(mols_targetsim_Morgen, key=lambda x: x[1], reverse=True)
result_Morgen = sorted_mols_targetsim_Morgen[:300]

#Morgen
print("### Morgen ###")
print()
n = 1
for m,fp in result_Morgen:
  for i in range(1,16):   
    if m in globals()["mol"+str(i)]:
     print(str(n) + "_" + globals()["NAME"+str(i)] + "_" + m.GetProp("Name") + ": " + str(fp))
    else:
     continue
  n = n + 1
print()
print('Completed!')
print()
print('Writing Outputfile...')
out_Morgen = []
for mol,sim in result_Morgen:
    out_Morgen.append(mol)
    mol.SetProp('Morgen_Similarity', str(sim))
writer = Chem.SDWriter('result_Morgen.sdf')
for i,m in enumerate(out_Morgen):
    writer.write(m)
writer.close()
print('Finished!')
