import os
import sys
import shutil

def generate_scan_files(path,path_ik_s, path_tera_opt):
    info = read_ik_s_file(path_ik_s)
    for i, inf in enumerate(info):
        try:
            ik = inf[0]
            path_ik = os.path.join(path, ik)
            if not os.path.exists(path_ik):
                os.mkdir(path_ik)
            else:
                path_ik = os.path.join(path, ik + '_' + str(i))
                if not os.path.exists(path_ik):
                    os.mkdir(path_ik)
            scan_term = inf[2]
            charge = inf[1]
            start_s = generate_scan(charge,scan_term)

            path_start = os.path.join(path_ik, 'scan.start')
            with open(path_start, 'w')as f:
                f.write(start_s)

            path_xyz_o = os.path.join(path_tera_opt, ik, 'scr.original', 'optim.xyz')
            xyz_s = get_scan_xyz(path_xyz_o)
            path_xyz = os.path.join(path_ik, 'scan.xyz')
            with open(path_xyz, 'w') as f:
                f.write(xyz_s)

        #path_xyz_dest = os.path.join(path_ik, 'scan.xyz')
        #shutil.copyfile(path_xyz, path_xyz_dest)
        #path_xyz_new = os.path.join(path_ik, 'xtbopt.xyz')
        #path_xyz_modi = os.path.join(path_ik, 'scan.xyz')
        #os.system('mv {} {}'.format(path_xyz_new, path_xyz_modi))
        except:
            print(ik, 'is error')

def get_scan_xyz(path):
    with open(path,'r')as f:
        lines = f.readlines()
    atom_num = int(lines[0].split()[0])+2
    xyz_lines = lines[-atom_num:]
    s = ''
    for line in xyz_lines:
        s +=line
    return s
def generate_scan(charge,scan_term):
    spin = 1
    method = 'b3lyp'
    if int(charge) % 2 != 0:
        spin = 1
        method = 'rob3lyp'
    s_t = scan_term.split('_')
    s_t = [str(int(i)+1) for i in s_t]
    s_t = '_'.join(s_t)
    s = '$multibasis\n'
    s += 'Br   6-311g*\n'
    s += '$end\n'
    s += 'coordinates scan.xyz\n'
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
    s += 'end\n\n'
    s += '$constraint_scan\n'
    s += 'dihedral 0.0, 360 13 {}\n'.format(s_t)
    s += '$end'
    return s

def read_ik_s_file(path_ik_s):
    info = []
    with open(path_ik_s,'r') as f:
        lines = f.readlines()

    for line in lines:
        inf = line.split(',')
        if len(inf)>0:
            info.append(inf)
    return info

if __name__ == "__main__":
    # the path of source
    path_tera_opt = '/nfs2/zcc/new_comparing_process/tera_opt'
    path_ik_s = '/nfs2/zcc/new_comparing_process/ik_smiles.stxm'
    path_tera_scan = '/nfs2/zcc/new_comparing_process/tera_scan'
    if not os.path.exists(path_tera_scan):
        os.mkdir(path_tera_scan)

    generate_scan_files(path_tera_scan, path_ik_s,path_tera_opt )


