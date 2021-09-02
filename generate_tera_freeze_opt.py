
import os
import sys
import shutil
from rdkit import Chem
def read_xtbscan_log(path):
    with open(path,'r') as f:
        lines = f.readlines()
    a_num = int(lines[0].split()[0])
    confs = []
    i = 0
    while i < len(lines):
        conf = lines[i:a_num+2+i]
        s = ''
        for cn in conf:
            s += cn
        confs.append(s)
        i += a_num+2
    return confs

def generate_start(charge,dih_l,n_xyz):  # int list xyz_name
    spin = 1
    method = 'b3lyp'
    if int(charge) % 2 != 0:
        spin = 1
        method = 'rob3lyp'
    s = '$multibasis\n'
    s += 'Br   6-311g*\n'
    s += '$end\n'
    s += 'coordinates {}\n'.format(n_xyz)
    s += 'charge {}\n'.format(charge)
    s += 'spinmult {}\n'.format(spin)
    s += 'min_maxiter 1000\n'
    s += 'run minimize\n'
    s += 'basis 6-31g**\n'
    s += 'maxit 500\n'
    s += 'convthre 1.0e-4\n'
    s += 'new_minimizer yes\n'
    s += 'method {}\n'.format(method)
    s += 'dftd no\n'
    s += 'end\n'
    s += '$constraint_freeze\n'
    for dih in dih_l:
        s += 'dihedral {}\n'.format(dih)  # 1_2_3_4
    s += '$end'
    return s

def read_ik_s_file(path_ik_s):
    info = {}
    with open(path_ik_s,'r') as f:
        lines = f.readlines()

    for line in lines:
        inf = line.split(',')
        if len(inf)>0:
            #print(inf)
            info[inf[0]] = [inf[1],inf[-1]]
    return info
def get_all_dihedral(smiles):
    dihedral = []
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
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
                scan = [str(head+1), str(center_atoms[0]+1), str(center_atoms[-1]+1), str(tail+1)]
                scan_str = '_'.join(scan)
                #dihe = '{},{},{},{}'.format(ik, charge, scan_str,smiles)
                dihedral.append(scan_str)
    return dihedral
def get_all_dihedral_all_s_d_f(smiles):
    dihedral = []
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    for bond in bonds:
        #print(str(bond.GetBondType())=='SINGLE',not bond.IsInRing())
        if str(bond.GetBondType()) == 'SINGLE':
            if not bond.IsInRing():
                #print('chain single')
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
                scan = [str(head+1), str(center_atoms[0]+1), str(center_atoms[-1]+1), str(tail+1)]
                #scan = [str(head), str(center_atoms[0] ), str(center_atoms[-1] ), str(tail )]
                scan_str = '_'.join(scan)
                #print(scan_str)
                #dihe = '{},{},{},{}'.format(ik, charge, scan_str,smiles)
                dihedral.append(scan_str)

            if bond.IsInRing():
                #print('ring single')
                atom_0 = bond.GetBeginAtomIdx()
                atom_1 = bond.GetEndAtomIdx()
                center_atoms = [atom_0, atom_1]
                nei_atoms_0 = atoms[center_atoms[0]].GetNeighbors()
                nei_atoms_0_idx = [nei_atom.GetIdx() for nei_atom in nei_atoms_0]
                if len(set(nei_atoms_0_idx))>1:
                    heads = list(set(nei_atoms_0_idx) - set(center_atoms))
                    head_0 = 0
                    for head in heads:
                        if atoms[head].IsInRing():
                            head_0 = head
                else:
                    continue
                nei_atoms_1 = atoms[center_atoms[-1]].GetNeighbors()
                nei_atoms_1_idx = [nei_atom.GetIdx() for nei_atom in nei_atoms_1]
                if len(set(nei_atoms_1_idx))>1:
                    tails = list(set(nei_atoms_1_idx) - set(center_atoms))
                    tail_0 = 0
                    for tail in tails:
                        if atoms[tail].IsInRing():
                            tail_0 = tail
                else:
                    continue
                scan = [str(head_0 + 1), str(center_atoms[0] + 1), str(center_atoms[-1] + 1), str(tail_0 + 1)]
                #scan = [str(head_0), str(center_atoms[0]), str(center_atoms[-1]), str(tail_0)]
                scan_str = '_'.join(scan)
                #print(scan_str)
                dihedral.append(scan_str)

    return dihedral
# generate terachem files from xtb scan
def generate_start_xyz_from_scan(path_scan,path_xtb_sp,smiles_file):
    base_info = read_ik_s_file(smiles_file)  # dict

    iks = os.listdir(path_scan)
    for ik in iks:
        path_xtbscan = os.path.join(path_scan,ik,'xtbscan.log')
        if os.path.exists(path_xtbscan):
            new_path_ik = os.path.join(path_xtb_sp,ik)
            if not os.path.exists(new_path_ik):
                os.mkdir(new_path_ik)
            #else:
                #new_path_ik = os.path.join(path_xtb_sp, ik+'_'+str(i))
            confs_xyz = read_xtbscan_log(path_xtbscan)
            for i,conf in enumerate(confs_xyz):
                n_xyz = 'strain_'+str(i)+'.xyz'
                n_start = 'strain_' + str(i) + '.start'
                path_xyz = os.path.join(new_path_ik,n_xyz)
                with open(path_xyz,'w') as f:
                    f.write(conf)
                path_start = os.path.join(new_path_ik,n_start)
                ik = ik.split("_")[0]
                charge = base_info[ik][0]
                smiles = base_info[ik][-1]
                all_dihe = get_all_dihedral_all_s_d_f(smiles)
                #print(all_dihe)
                s_start = generate_start(charge,all_dihe,n_xyz)
                with open(path_start,'w') as f:
                    f.write(s_start)

def generate_start_xyz_from_opt(path_scan,path_xtb_sp,smiles_file):
    base_info = read_ik_s_file(smiles_file)

    iks = os.listdir(path_scan)
    for ik in iks:
        path_xtbscan = os.path.join(path_scan,ik,'xtbopt.xyz')
        if os.path.exists(path_xtbscan):
            new_path_ik = os.path.join(path_xtb_sp,ik)
            if not os.path.exists(new_path_ik):
                os.mkdir(new_path_ik)
            #else:
                #new_path_ik = os.path.join(path_xtb_sp, ik+'_'+str(i))
            confs_xyz = read_xtbscan_log(path_xtbscan)
            for i,conf in enumerate(confs_xyz):
                try:
                    n_xyz = 'strain_' + str(i) + '.xyz'
                    n_start = 'strain_' + str(i) + '.start'
                    path_xyz = os.path.join(new_path_ik, n_xyz)
                    with open(path_xyz, 'w') as f:
                        f.write(conf)
                    path_start = os.path.join(new_path_ik, n_start)
                    ik = ik.split("_")[0]
                    charge = base_info[ik][0]
                    smiles = base_info[ik][-1]
                    all_dihe = get_all_dihedral_all_s_d_f(smiles)

                    s_start = generate_start(charge, all_dihe, n_xyz)
                    with open(path_start, 'w') as f:
                        f.write(s_start)
                except:
                    print(ik, ' is error')


if __name__ == "__main__":
    '''
    smiles_file = '/nfs2/zcc/new_comparing_process/BrccCl/ik_smiles.stxm'

    path_xtb_scan = '/nfs2/zcc/new_comparing_process/BrccCl/xtb_scan'
    #path_xtb_scan = '/nfs2/zcc/new_comparing_process/results/lines_moles_result/xtb_opt'
    path_xtb_tera_starin = '/nfs2/zcc/new_comparing_process/BrccCl/xtb_scan_tera_strain'
    if not os.path.exists(path_xtb_tera_starin):
        os.mkdir(path_xtb_tera_starin)

    generate_start_xyz_from_scan(path_xtb_scan, path_xtb_tera_starin,smiles_file)
    '''
    #'''
    smiles_file = '/nfs2/zcc/new_comparing_process/frag_ik_smiles.stxm'

    path_xtb_opt = '/nfs2/zcc/new_comparing_process/frag_xtb_opt'
    # path_xtb_scan = '/nfs2/zcc/new_comparing_process/results/lines_moles_result/xtb_opt'
    path_xtb_tera_starin = '/nfs2/zcc/new_comparing_process/frag_xtb_opt_tera_strain'
    if not os.path.exists(path_xtb_tera_starin):
        os.mkdir(path_xtb_tera_starin)

        generate_start_xyz_from_opt(path_xtb_opt, path_xtb_tera_starin, smiles_file)
    #'''