#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Aug  9 14:25:59 2021

@author: silicon
"""


from rdkit import Chem
from rdkit.Chem import AllChem
import os


def get_mols_info(mol_smiles,fun_ground,patt_smiles,c_atoms):
    
    info = []
    mol = Chem.MolFromSmiles(mol_smiles)
    mol_atoms = mol.GetAtoms()
    mol_atoms[c_atoms[0]].SetProp('center_atom','true')
    mol_atoms[c_atoms[1]].SetProp('center_atom','true')
    
    for fg in fun_ground:
        try:
            repl = Chem.MolFromSmiles(fg)
            patt = Chem.MolFromSmarts(patt_smiles)
            new_mol = AllChem.ReplaceSubstructs(mol, patt, repl, replacementConnectionPoint=0)[0]
            atoms = new_mol.GetAtoms()
            center_atoms = []
            for i, at in enumerate(atoms):
                if at.HasProp('center_atom'):
                    center_atoms.append(i)
            nei_atoms_0 = atoms[center_atoms[0]].GetNeighbors()
            nei_atoms_0_idx = [nei_atom.GetIdx() for nei_atom in nei_atoms_0]
            head = list(set(nei_atoms_0_idx)-set(center_atoms))[0]
            
            nei_atoms_1 = atoms[center_atoms[-1]].GetNeighbors()
            nei_atoms_1_idx = [nei_atom.GetIdx() for nei_atom in nei_atoms_1]
            tail = list(set(nei_atoms_1_idx)-set(center_atoms))[0]
            scan = [str(head),str(center_atoms[0]),str(center_atoms[-1]),str(tail)]
            new_smiles = Chem.MolToSmiles(new_mol)
            info.append([scan,new_smiles])
        except:
            print(fg, 'is error')
            continue
    return info

def write_xtb_opt_job(ik, mol, path):
    if not os.path.exists(path):
        os.mkdir(path)
    fd = os.path.join(path,ik)
    if not os.path.exists(fd):
        os.mkdir(fd)
    filename = os.path.join(fd,'original.xyz')
    Chem.MolToXYZFile(mol, filename)

def write_tera_opt_job(ik, mol, path):
    if not os.path.exists(path):
        os.mkdir(path)
    fd = os.path.join(path,ik)
    if not os.path.exists(fd):
        os.mkdir(fd)
    filename = os.path.join(fd,'original.xyz')
    Chem.MolToXYZFile(mol, filename)
    start_file = os.path.join(fd,'opt.start')
    atoms = mol.GetAtoms()
    charge = 0

    for atom in atoms:
        charge += atom.GetFormalCharge()
    spin = 1
    method = 'b3lyp'
    if charge%2 != 0:
        spin = 1
        method = 'rob3lyp'
    s = '$multibasis\n'
    s+='Br   6-311g*\n'
    s+='$end\n'
    s+='coordinates original.xyz\n'
    s+='charge {}\n'.format(charge)
    s+='spinmult {}\n'.format(spin)
    s+='min_maxiter 1000\n'
    s+='run minimize\n'
    s+='basis 6-31g**\n'
    s+='maxit 500\n'
    s+='convthre 1.0e-4\n'
    s+='new_minimizer yes\n'
    s+='method {}\n'.format(method)
    s+='dftd no\n'
    s+='end\n'
    
    with open(start_file,'w') as f:
        f.write(s)
#
def generate_file(mol_smiles,fun_ground, patt_smiles,c_atoms, smiles_file, files_path):
    info =get_mols_info(mol_smiles,fun_ground, patt_smiles,c_atoms) 
     
    s = ''
    for ff in info:
        scan = ff[0] 
        scan_str = '_'.join(scan)
        smiles = ff[1].split('\n')[0]
        
        mol = Chem.MolFromSmiles(smiles) 
        ik = Chem.MolToInchiKey(mol)
        atoms = mol.GetAtoms()
        charge = 0
        for atom in atoms:
            charge += atom.GetFormalCharge()
        s += '{},{},{},{}\n'.format(ik,charge,scan_str,smiles)
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        
        xtb_sp_path = os.path.join(files_path,'xtb_opt')
        if not os.path.exists(xtb_sp_path):
            os.mkdir(xtb_sp_path)
        write_xtb_opt_job(ik, mol, xtb_sp_path)

        tera_sp_path = os.path.join(files_path,'tera_opt')
        if not os.path.exists(tera_sp_path):
            os.mkdir(tera_sp_path)
        write_tera_opt_job(ik, mol, tera_sp_path)
        
    with open(smiles_file,'w') as f:
        f.write(s)
def  generate_file_1(smiles,out_files_path,smiles_file):
    smiles = smiles.split('\n')[0]
    mol = Chem.MolFromSmiles(smiles)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()

    #mol = Chem.MolFromSmiles(smiles)
    ik = Chem.MolToInchiKey(mol)
    charge = 0
    for atom in atoms:
        charge += atom.GetFormalCharge()
    s = ''
    for bond in bonds:
        #print(str(bond.GetBondType())=='SINGLE',not bond.IsInRing())
        if str(bond.GetBondType()) == 'SINGLE':
            if not bond.IsInRing():
                atom_0 = bond.GetBeginAtomIdx()
                atom_1 = bond.GetEndAtomIdx()
                center_atoms = [atom_0,atom_1]
                nei_atoms_0 = atoms[center_atoms[0]].GetNeighbors()
                nei_atoms_0_idx = [nei_atom.GetIdx() for nei_atom in nei_atoms_0]
                if len(set(nei_atoms_0_idx))>1:
                    head = list(set(nei_atoms_0_idx) - set(center_atoms))[0]
                else:
                    continue
                nei_atoms_1 = atoms[center_atoms[-1]].GetNeighbors()
                nei_atoms_1_idx = [nei_atom.GetIdx() for nei_atom in nei_atoms_1]

                if len(set(nei_atoms_1_idx))>1:
                    tail = list(set(nei_atoms_1_idx) - set(center_atoms))[0]
                else:
                    continue
                scan = [str(head), str(center_atoms[0]), str(center_atoms[-1]), str(tail)]
                scan_str = '_'.join(scan)
                s += '{},{},{},{}\n'.format(ik, charge, scan_str,smiles)
                #s += '{},{},{},{}'.format(ik, charge, scan_str, smiles)
                #print('i am bad')
                #print(s)
    print(s)
    with open(smiles_file,'a') as f:
        f.write(s)

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    xtb_sp_path = os.path.join(out_files_path, 'frag_xtb_opt')
    if not os.path.exists(xtb_sp_path):
        os.mkdir(xtb_sp_path)
    write_xtb_opt_job(ik, mol, xtb_sp_path)

    tera_sp_path = os.path.join(out_files_path, 'frag_tera_opt')
    if not os.path.exists(tera_sp_path):
        os.mkdir(tera_sp_path)
    write_tera_opt_job(ik, mol, tera_sp_path)

def read_smiles(file_path):
    info = []  # [[ik.smiles]]
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        inf = line.split(',')
        if len(inf) > 0:
            info.append(inf)
    return info
if __name__ == "__main__":

    # the mothod 1
    '''
    mol_smiles = 'BrCCBr'
    #mol_smiles = 'c1ccccc1CC(N)C(=O)O'
    fun_ground = ['Br','CC','CCC','CCCC','C(=O)O','C[O-]','CC(CC)C=O','C1CCCCC1Br','CC(O)CCC','Cc1cscc1',]
    patt_smiles = 'Br'   # 需要替换官能团
    c_atoms = [1,2]
    # the path of generating ik_smiles.stxm
    smiles_file = '/nfs2/zcc/new_comparing_process/ik_smiles.stxm'
    # the path of putting opt files
    out_files_path = '/nfs2/zcc/new_comparing_process'
    generate_file(mol_smiles,fun_ground, patt_smiles, c_atoms, smiles_file, out_files_path)
    '''
    # THE METHOD 2
    #'''
    #col_smiles = ['COC(=O)Cc1cccc(C[C@H](N)C(=O)O)c1','FC(F)(F)c1cccc(-c2nn[nH]n2)c1','CCC(N)C(=O)N1CCCC1']
    path_ik_smiles = '/nfs2/zcc/new_comparing_process/frag.stxm'
    ik_smiles = read_smiles(path_ik_smiles)
    #for smiles in col_smiles:
    for ik_s in ik_smiles[0:1000]:
        #smiles = 'OCC(O)C(O)C(O)C(O)C=O'
        smiles = ik_s[-1]
        out_files_path = '/nfs2/zcc/new_comparing_process'
        smiles_file = '/nfs2/zcc/new_comparing_process/frag_ik_smiles.stxm'
        generate_file_1(smiles, out_files_path, smiles_file)
    #'''



    






       
