import os
import sys
import shutil

def generate_scan_files(path,path_ik_s, path_xtb_opt):
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
            inp_s = generate_inp(scan_term)
            path_inp = os.path.join(path_ik, 'scan.inp')
            with open(path_inp, 'w')as f:
                f.write(inp_s)
            path_xyz = os.path.join(path_xtb_opt, ik, 'xtbopt.xyz')
            path_xyz_dest = os.path.join(path_ik, 'scan.xyz')
            shutil.copyfile(path_xyz, path_xyz_dest)

        #path_xyz_new = os.path.join(path_ik, 'xtbopt.xyz')
        #path_xyz_modi = os.path.join(path_ik, 'scan.xyz')
        #os.system('mv {} {}'.format(path_xyz_new, path_xyz_modi))
        except:
            print(ik, 'is error')

def generate_inp(scan_term):
    s_t = scan_term.split('_')
    s_t = [str(int(i) + 1) for i in s_t]
    s = '$constrain\n'
    s += 'force constant=1.0\n'
    s += 'dihedral: {},{},{},{}, 0.0\n'.format(s_t[0],s_t[1],s_t[2],s_t[3])
    s += '$scan\n'
    s += '1: 0.0, 360, 13\n'
    s += '$end\n'
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
    path_xtb_opt = '/nfs2/zcc/new_comparing_process/frag_xtb_opt'
    path_ik_s = '/nfs2/zcc/new_comparing_process/frag_ik_smiles.stxm'
    path_xtb_scan = '/nfs2/zcc/new_comparing_process/frag_xtb_scan'
    if not os.path.exists(path_xtb_scan):
        os.mkdir(path_xtb_scan)

    generate_scan_files(path_xtb_scan, path_ik_s, path_xtb_opt)


