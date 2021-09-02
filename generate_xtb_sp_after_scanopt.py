import os
import sys
import shutil

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

def generate_start(charge,n_xyz):
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
    s += 'run energy\n'
    s += 'basis 6-31g**\n'
    s += 'convthre 1.0e-4\n'
    s += 'method {}\n'.format(method)
    s += 'dftd no\n'
    s += 'end\n'
    return s

def read_ik_s_file(path_ik_s):
    info = {}
    with open(path_ik_s,'r') as f:
        lines = f.readlines()

    for line in lines:
        inf = line.split(',')
        if len(inf)>0:
            info[inf[0]] = inf[1]
    return info

# generate terachem files from xtb scan
def generate_start_xyz(path_scan,path_xtb_sp,smiles_file):
    base_info = read_ik_s_file(smiles_file)

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
                n_xyz = 'sp_'+str(i)+'.xyz'
                n_start = 'sp_' + str(i) + '.start'
                path_xyz = os.path.join(new_path_ik,n_xyz)
                with open(path_xyz,'w') as f:
                    f.write(conf)
                path_start = os.path.join(new_path_ik,n_start)
                ik = ik.split("_")[0]
                charge = base_info[ik]
                s_start = generate_start(charge,n_xyz)
                with open(path_start,'w') as f:
                    f.write(s_start)
# generate terachem files from xtb opt
def generate_start_xyz_1(path_scan,path_xtb_sp,smiles_file):
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
                    n_xyz = 'sp_' + str(i) + '.xyz'
                    n_start = 'sp_' + str(i) + '.start'
                    path_xyz = os.path.join(new_path_ik, n_xyz)
                    with open(path_xyz, 'w') as f:
                        f.write(conf)
                    path_start = os.path.join(new_path_ik, n_start)
                    ik = ik.split("_")[0]

                    charge = base_info[ik]

                    s_start = generate_start(charge, n_xyz)
                    with open(path_start, 'w') as f:
                        f.write(s_start)
                except:
                    print(ik, ' is error')
                    continue


if __name__ == "__main__":
    #generate terachem files from xtb scan
    '''
    smiles_file = '/nfs2/zcc/new_comparing_process/ik_smiles.stxm'

    path_xtb_scan = '/nfs2/zcc/new_comparing_process/xtb_scan'
    #path_xtb_scan = '/nfs2/zcc/new_comparing_process/results/lines_moles_result/xtb_opt'
    path_xtb_tera_sp = '/nfs2/zcc/new_comparing_process/xtb_opt_tera_sp'
    if not os.path.exists(path_xtb_tera_sp):
        os.mkdir(path_xtb_tera_sp)

    generate_start_xyz(path_xtb_scan, path_xtb_tera_sp,smiles_file)
    '''

    # use for terachem sp files based on xtb opt files
    #'''
    smiles_file = '/nfs2/zcc/new_comparing_process/frag_ik_smiles.stxm'

    path_xtb_opt = '/nfs2/zcc/new_comparing_process/frag_xtb_opt'
    # path_xtb_scan = '/nfs2/zcc/new_comparing_process/results/lines_moles_result/xtb_opt'
    path_xtb_tera_sp = '/nfs2/zcc/new_comparing_process/frag_xtb_opt_tera_sp'
    if not os.path.exists(path_xtb_tera_sp):
        os.mkdir(path_xtb_tera_sp)

    generate_start_xyz_1(path_xtb_opt, path_xtb_tera_sp, smiles_file)
    #'''
