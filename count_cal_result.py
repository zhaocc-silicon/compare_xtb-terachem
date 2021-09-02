import os
import sys
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import numpy as np

def get_xtb_energy(xtb_sp_path):
    iks = os.listdir(xtb_sp_path)
    xtb_info = {}
    for ik in iks:
        try:
            en_c = []
            path_ik = os.path.join(xtb_sp_path, ik)
            fs = os.listdir(path_ik)
            all_scr = []
            for ff in fs:
                if 'scr.' in ff:
                    all_scr.append(ff)
            #new_scr = sorted(all_scr, key=lambda i: i[-1])
            new_scr = sorted(all_scr, key=lambda i: int(i.split("_")[-1]))
            for scr in new_scr:
                path_xyz = os.path.join(path_ik, scr, 'xyz.xyz')
                with open(path_xyz, 'r') as f:
                    lines = f.readlines()
                energy = lines[1].split()[0]
                path_mol = os.path.join(path_ik, scr, 'xyz.mol')
                os.system('obabel -ixyz {} -omol -O {}'.format(path_xyz, path_mol))
                mol_i = Chem.MolFromMolFile(path_mol)
                en_c.append([energy, mol_i])
            xtb_info[ik] = en_c
        except:
            print(ik, 'xtb sp calculated error')
    return xtb_info  #{ik:[[energy,moli],...],..}

def get_xtb_energy_1(xtb_sp_path):
    iks = os.listdir(xtb_sp_path)
    xtb_info = {}
    for ik in iks:
        try:
            en_c = {}
            path_ik = os.path.join(xtb_sp_path, ik)
            fs = os.listdir(path_ik)
            all_scr = []
            for ff in fs:
                if 'scr.' in ff:
                    all_scr.append(ff)
            #new_scr = sorted(all_scr, key=lambda i: i[-1])
            new_scr = sorted(all_scr, key=lambda i: int(i.split("_")[-1]))
            for i, scr in enumerate(new_scr):
                path_xyz = os.path.join(path_ik, scr, 'xyz.xyz')
                with open(path_xyz, 'r') as f:
                    lines = f.readlines()
                energy = lines[1].split()[0]
                path_mol = os.path.join(path_ik, scr, 'xyz.mol')
                os.system('obabel -ixyz {} -omol -O {}'.format(path_xyz, path_mol))
                mol_i = Chem.MolFromMolFile(path_mol)
                en_c[i] = ([energy, mol_i])
            xtb_info[ik] = en_c
        except:
            print(ik, 'xtb sp calculated error')
    return xtb_info  #{ik:[[energy,moli],...],..}

def get_tera_constrain_energy(xtb_scan_tera_strain):
    iks = os.listdir(xtb_scan_tera_strain)
    tera_info = {}
    for ik in iks:
        #try:
            en_c = {}
            path_ik = os.path.join(xtb_scan_tera_strain, ik)
            fs = os.listdir(path_ik)
            all_scr = []
            for ff in fs:
                if 'scr.' in ff:
                    all_scr.append(ff)

            #new_scr = sorted(all_scr, key=lambda i: i[-1])
            new_scr = sorted(all_scr, key=lambda i: int(i.split("_")[-1]))
            for ii, scr in enumerate(new_scr):
                path_xyz = os.path.join(path_ik, scr, 'optim.xyz')
                #with open(path_xyz, 'r') as f:
                    #lines = f.readlines()
                #en_conf = read_tera_opt_xyz(path_xyz)
                try:
                    en_confs = read_tera_opt_xyz(path_xyz)
                except:
                    print(path_xyz, '  error')
                    continue
                for i, en_f in enumerate(en_confs):
                    energy = en_f[0]
                    xyz_s = en_f[1]
                    path_xyz = os.path.join(path_ik, scr, str(i) + '.xyz')
                    with open(path_xyz, 'w') as f:
                        f.write(xyz_s)
                    path_mol = os.path.join(path_ik, scr, str(i) + '.mol')
                    os.system('obabel -ixyz {} -omol -O {}'.format(path_xyz, path_mol))
                    mol_i = Chem.MolFromMolFile(path_mol)
                    en_c[ii] = ([energy, mol_i])
            tera_info[ik] = en_c    # {ik:[[energy,mol]]}
    return tera_info


def read_tera_scan_xyz(path):
    with open(path,'r') as f:
        lines = f.readlines()
    a_num = int(lines[0].split()[0])
    en_confs = []
    i = 0
    while i < len(lines):
        conf = lines[i:a_num+2+i]
        s = ''
        for cn in conf:
            s += cn
            if 'Energy' in cn:
                energy = float(cn.split()[4])
        en_confs.append([energy,s])
        i += a_num+2
    return en_confs
# get info from tera_scan
def get_tera_mol(tara_scan_path):
    iks = os.listdir(tara_scan_path)
    tera_info = {}
    for ik in iks:
        try:
            path_ik = os.path.join(tara_scan_path, ik)
            fs = os.listdir(path_ik)
            en_c = []
            for ff in fs:
                if 'scr.' in ff:
                    path_scan_opt = os.path.join(path_ik, ff, 'scan_optim.xyz')
                    if os.path.exists(path_scan_opt):
                        #print(path_scan_opt)
                        try:
                            en_confs = read_tera_scan_xyz(path_scan_opt)
                        except:
                            print(ik, 'none scan_optim.xyz')
                            continue
                        for i, en_f in enumerate(en_confs):
                            energy = en_f[0]
                            xyz_s = en_f[1]
                            path_xyz = os.path.join(path_ik, ff, str(i) + '.xyz')
                            with open(path_xyz, 'w') as f:
                                f.write(xyz_s)
                            path_mol = os.path.join(path_ik, ff, str(i) + '.mol')
                            os.system('obabel -ixyz {} -omol -O {}'.format(path_xyz, path_mol))
                            mol_i = Chem.MolFromMolFile(path_mol)
                            en_c.append([energy, mol_i])
                        break
            tera_info[ik] = en_c
        except:
            print(ik,'tera scan calculated error')
    return tera_info  #{ik:[[energy,moli],...],..}

# plt de_i, rmse
def plot_picture(path_plt,diff_energies,conf_rmse):

    fig = plt.figure(figsize=(15,9))
    length = len(diff_energies)
    step = int(360/(length-1))
    #print(step)
    x = list(range(0,360+step,step))
    ax1 = fig.add_subplot(111)
    y1 = np.array(diff_energies)*627.5094706
    y2 = conf_rmse
    lns1 = ax1.plot(x,y1,'*:r',label='energy difference')
    ax1.set_ylabel('diff_energy(kcal/mol)',fontsize=16)
    ax1.set_xlabel('dihedral',fontsize=16)

    ax2 = ax1.twinx()
    lns2 = ax2.plot(x,y2,'-o',color='b',label='RMSE')
    ax2.set_ylabel('rmse(A)',fontsize=16)
    ax2.set_xlabel('dihedral',fontsize=16)
    #labels = x
    lns = lns1+lns2
    labels = [l.get_label() for l in lns]
    ax1.legend(lns, labels,loc=0,fontsize=16)
    plt.xticks(x,fontsize=16)
    fig_path = os.path.join(path_plt,'en_rms.png')
    plt.savefig(fig_path)
    #plt.show()
    plt.close()


# get info from tera_scan
def get_tera_from_opt(tara_opt_path):
    iks = os.listdir(tara_opt_path)
    tera_info = {}
    for ik in iks:
        try:
            path_ik = os.path.join(tara_opt_path, ik)
            fs = os.listdir(path_ik)
            en_c = []
            for ff in fs:
                if 'scr.' in ff:
                    path_opt = os.path.join(path_ik, ff, 'optim.xyz')
                    if os.path.exists(path_opt):
                        #print(path_scan_opt)
                        try:
                            en_conf = read_tera_opt_xyz(path_opt)
                        except:
                            print(ik, 'none optim.xyz')
                            continue
                        for i, en_f in enumerate(en_conf):
                            energy = en_f[0]
                            xyz_s = en_f[1]
                            path_xyz = os.path.join(path_ik, ff, str(i) + '.xyz')
                            with open(path_xyz, 'w') as f:
                                f.write(xyz_s)
                            path_mol = os.path.join(path_ik, ff, str(i) + '.mol')
                            os.system('obabel -ixyz {} -omol -O {}'.format(path_xyz, path_mol))
                            mol_i = Chem.MolFromMolFile(path_mol)
                            en_c.append([energy, mol_i])
                        break
            tera_info[ik] = en_c
        except:
            print(ik,'tera scan calculated error')
    return tera_info  #{ik:[[energy,moli],...],..}

def read_tera_opt_xyz(path):
    info = []
    with open(path, 'r') as f:
        lines = f.readlines()
    a_num = int(lines[0].split()[0])
    energy = float(lines[-(a_num+1)].split()[0])
    s = ''
    for line in lines[-(a_num+2):]:
        s += line
    info.append([energy,s])
    #print(energy)
    #print(s)
    return info

# plt PES, rmse
def plot_picture_1(path_plt,energies,conf_rmse):
    energies = np.array(energies)
    xtb_energy = energies[:,0]
    xtb_min = np.min(xtb_energy)
    xtb_energy =xtb_energy - xtb_min

    tera_energy = energies[:,1]
    tera_min = np.min(tera_energy)
    tera_energy = tera_energy - tera_min

    fig = plt.figure(figsize=(15,9))
    length = len(energies)
    step = int(360/(length-1))
    #print(step)
    x = list(range(0,360+step,step))
    ax1 = fig.add_subplot(111)
    y1 = xtb_energy*627.5094706
    y3 = tera_energy*627.5094706
    y2 = conf_rmse
    #lns1 = ax1.plot(x,y1,'*:r',label='PES of xtb')
    lns1 = ax1.plot(x, y1, '-o',color='r',label='PES of xtb')
    #lns3 = ax1.plot(x,y3,'.:r',label='PES of tera')
    lns3 = ax1.plot(x, y3, '-o',color='m', label='PES of tera')
    ax1.set_ylabel('energy(kcal/mol)',color='r',fontsize=16)
    ax1.set_xlabel('dihedral(degree)',color='k',fontsize=16)

    ax2 = ax1.twinx()
    lns2 = ax2.plot(x,y2,'--^',color='b',label='RMSE')
    ax2.set_ylabel('rmse(A)',color='b',fontsize=16)
    ax2.set_xlabel('dihedral',fontsize=16)
    #labels = x
    lns = lns1+lns3+lns2
    labels = [l.get_label() for l in lns]
    ax1.legend(lns, labels,loc=0,fontsize=16)
    plt.xticks(x,fontsize=16)
    fig_path = os.path.join(path_plt,'en_rms.png')
    plt.savefig(fig_path)
    #plt.show()
    plt.close()

# plt sctter diagram
def plot_picture_2(path_plt,all_energies):
    #print("i want to plot sctter")
    colors = ['k','g','b','c','y','m','r','pink','chocolate']
    plt.figure(figsize=(8, 6))
    ax = plt.subplot()
    ax.plot((0,1),(0,1),transform=ax.transAxes, ls='--',c='k')
    i = 0
    print(all_energies.values())
    en_values = all_energies.values()
    values = np.array(list(en_values))
    #max_en = np.max(values)
    #min_en = np.min(values)
    #print('max_en',max_en)
    #print('min_en',min_en)
    #space = int((max_en - min_en)*627.5094706*1.2)
    #step = int(space/10)
    #space = step*10
    #ticks = list(range(0,space,step))
    ticks = list(range(0, 15, 3))
    for ik, energies in all_energies.items():
        energies = np.array(energies)
        #print(energies)
        xtb_energy = energies[:,0]*627.5094706
        xtb_min = np.min(xtb_energy)
        xtb_energy = xtb_energy - xtb_min

        tera_energy = energies[:,1]*627.5094706
        tera_min = np.min(tera_energy)
        tera_energy = tera_energy - tera_min

        # energies = np.array(energies)
        # xtb_energy = energies[:,0]
        # tera_energy = energies[:,1]
        if i < len(colors)-1:
            i+=1
        else:
            i = 0
        color_i = colors[i]
        #plt.scatter(tera_energy, xtb_energy,color=color_i, s=20, label='{}'.format(ik))
        plt.scatter(tera_energy, xtb_energy, s=20, color=color_i)
    plt.xlabel('b3lyp/6-31g** reference(kcal/mol)',fontsize=20)
    plt.ylabel('GFN2-xtb(kcal/mol)',fontsize=20)
    plt.xticks(ticks,size=20)
    plt.yticks(ticks,size=20)
    fig_path = os.path.join(path_plt, 'xtb_tera_scatter.png')
    plt.legend()
    plt.savefig(fig_path)
    # plt.show()
    plt.close()
# it cannot be plotted
def plot_opt_scatter(path_plt,all_energies):
    #print("i want to plot sctter")
    colors = ['k','g','b','c','y','m','r','pink','chocolate']
    plt.figure(figsize=(8, 6))
    ax = plt.subplot()
    ax.plot((0,1),(0,1),transform=ax.transAxes, ls='--',c='k')
    i = 0
    #print(all_energies.values())
    en_values = all_energies.values()
    values = np.array(list(en_values))
    max_en = np.max(values)
    min_en = np.min(values)
    print('max_en',max_en)
    print('min_en',min_en)
    space = int((max_en - min_en)*627.5094706*1.2)
    step = int(space/10)
    space = step*10
    ticks = list(range(0,space,step))
    #ticks = list(range(0, 15, 3))
    for ik, energies in all_energies.items():
        energies = np.array(energies)
        #print(energies)
        xtb_energy = energies[:,0]*627.5094706
        #xtb_min = np.min(xtb_energy)
        #xtb_energy = xtb_energy - xtb_min

        tera_energy = energies[:,1]*627.5094706
        #tera_min = np.min(tera_energy)
        #tera_energy = tera_energy - tera_min

        # energies = np.array(energies)
        # xtb_energy = energies[:,0]
        # tera_energy = energies[:,1]
        if i < len(colors)-1:
            i+=1
        else:
            i = 0
        color_i = colors[i]
        #plt.scatter(tera_energy, xtb_energy,color=color_i, label='{}'.format(ik))
        plt.scatter(tera_energy, xtb_energy, color=color_i)
    plt.xlabel('b3lyp/6-31g** reference(kcal/mol)',fontsize=20)
    plt.ylabel('GFN2-xtb(kcal/mol)',fontsize=20)
    plt.xticks(ticks,size=20)
    plt.yticks(ticks,size=20)
    fig_path = os.path.join(path_plt, 'xtb_tera_scatter.png')
    plt.legend()
    plt.savefig(fig_path)
    # plt.show()
    plt.close()
def plt_distri_energy(path_plt,energy):
    # energy --> list: [de0,de1,...]
    #plt.hist(energy, bins=15, range=(-15, 15), density=True)
    plt.figure(figsize=(8, 6))
    plt.xlabel('dE(kcal/mol)', fontsize=20)
    plt.xticks(size=20)
    plt.ylabel('Density', fontsize=20)
    plt.yticks(size=20)
    plt.hist(energy, bins=20, range=(-10, 10), density=True,color='r' ,histtype='step') # possiblity = 1*y
    path_energy = os.path.join(path_plt,'energy.png')
    plt.savefig(path_energy, dpi=300)
    plt.close()

def plt_distri_rmse(path_plt,rmse):
    # rmse --> list: [rm0,rm1,...]
    #plt.hist(energy, bins=15, range=(-15, 15), density=True)
    plt.figure(figsize=(8, 6))
    plt.xlabel('RMSE(A)', fontsize=20)
    plt.xticks(size=20)
    plt.ylabel('Density', fontsize=20)
    plt.yticks(size=20)
    plt.hist(rmse, bins=20, range=(0, 2.0), density=True, color='r', histtype='step')  # possiblity = 0.1*y
    path_energy = os.path.join(path_plt,'rmse.png')
    plt.savefig(path_energy, dpi=300)
    plt.close()

def write_basic_info(path_plt, energies, conf_rmse):
    s = 'xtb_energy(a.u.),tera_energy(a.u.),rmse(A)\n'
    for i, en in enumerate(energies):
        s+= '{},{},{}\n'.format(str(en[0]),str(en[1]),str(conf_rmse[i]))
    path_file = os.path.join(path_plt,'base_info.log')
    with open(path_file,'w') as f:
        f.write(s)

def compare_rmsd_energy(result_path,tara_scan_path,xtb_sp_path):   # the main funtion
    #info = {}
    all_energies = {}
    tera_infos = get_tera_mol(tara_scan_path)
    xtb_infos = get_xtb_energy(xtb_sp_path)
    diff_energies = []
    all_rmse = []
    for ik, xtb_value in xtb_infos.items():
        #print('start counting',ik,':')
        try:
            #diff_energies = []
            energies = []
            conf_rmse = []
            tera_value = tera_infos[ik]
            for i in range(len(xtb_value)):
                en_xtb = xtb_value[i][0]
                m_xtb = xtb_value[i][1]
                en_tera = tera_value[i][0]
                m_tera = tera_value[i][1]
                #print(Chem.MolToSmiles(m_xtb), Chem.MolToInchiKey(m_xtb))
                #print(Chem.MolToSmiles(m_tera), Chem.MolToInchiKey(m_tera))
                try:
                    rmse = AllChem.GetBestRMS(m_xtb, m_tera)
                except:
                    print('xtb:',Chem.MolToSmiles(m_xtb), Chem.MolToInchiKey(m_xtb))
                    print('tera:',Chem.MolToSmiles(m_tera), Chem.MolToInchiKey(m_tera))
                    print(ik, i,'conformations change too much')
                    break

                diff_en = (float(en_xtb) - float(en_tera)) * 627.5094706
                diff_energies.append(diff_en)
                conf_rmse.append(rmse)
                all_rmse.append(rmse)
                energies.append([float(en_xtb), float(en_tera)])

                #diff_en = float(en_xtb) - float(en_tera)
                #diff_energies.append(diff_en)
                #conf_rmse.append(rmse)

               # energies.append([float(en_xtb), float(en_tera)])
            if len(energies) >0:
                all_energies[ik] = energies
            path_plt = os.path.join(result_path, ik)
            if not os.path.exists((path_plt)):
                os.mkdir(path_plt)
            #plot_picture(path_plt, diff_energies, conf_rmse)
            plot_picture_1(path_plt, energies, conf_rmse)
            write_basic_info(path_plt, energies, conf_rmse)

        except:
            print(ik, 'is error')
        #info[ik] = [energies,conf_rmse]
    #'''
    #plt scatter
    path_plt = result_path
    #print(all_energies)
    plot_picture_2(path_plt, all_energies)
    plt_distri_energy(path_plt, diff_energies)
    plt_distri_rmse(path_plt, all_rmse)
    #'''
# compare opt
def compare_opt_energy_difference(result_path,tara_opt_path,xtb_sp_path):   # the main funtion
    #info = {}
    all_energies = {}
    tera_infos = get_tera_from_opt(tara_opt_path)
    #print(tera_infos)
    xtb_infos = get_xtb_energy(xtb_sp_path)
    #exit()
    diff_energies = []
    all_rmse = []
    #ii = 0   # if you want to restrict the molecules number, open it
    for ik, xtb_value in xtb_infos.items():
        #print('start counting',ik,':')
        #ii = ii+1
        #print(ii)
        try:
            #diff_energies = []
            energies = []
            conf_rmse = []
            tera_value = tera_infos[ik]
            for i in range(len(xtb_value)):
                en_xtb = xtb_value[i][0]
                m_xtb = xtb_value[i][1]
                en_tera = tera_value[i][0]
                m_tera = tera_value[i][1]
                #print(Chem.MolToSmiles(m_xtb), Chem.MolToInchiKey(m_xtb))
                #print(Chem.MolToSmiles(m_tera), Chem.MolToInchiKey(m_tera))
                try:
                    rmse = AllChem.GetBestRMS(m_xtb, m_tera)
                except:
                    print('xtb:',Chem.MolToSmiles(m_xtb), Chem.MolToInchiKey(m_xtb))
                    print('tera:',Chem.MolToSmiles(m_tera), Chem.MolToInchiKey(m_tera))
                    print(ik, i,'conformations change too much')
                    break
                diff_en = (float(en_xtb) - float(en_tera))*627.5094706
                diff_energies.append(diff_en)
                conf_rmse.append(rmse)
                all_rmse.append(rmse)
                energies.append([float(en_xtb), float(en_tera)])
            if len(energies) >0:
                all_energies[ik] = energies
            path_plt = os.path.join(result_path, ik)
            if not os.path.exists((path_plt)):
                os.mkdir(path_plt)
            #plot_picture(path_plt, diff_energies, conf_rmse)
            #plot_picture_1(path_plt, energies, conf_rmse)
            write_basic_info(path_plt, energies, conf_rmse)

        except:
            print(ik, 'is error')
        #if ii >=200:
            #print(ii)
            #break
        #info[ik] = [energies,conf_rmse]
    #'''
    #plt scatter
    path_plt = result_path
    #print(all_energies)
    #plot_opt_scatter(path_plt, all_energies)
    #'''
    plt_distri_energy(path_plt, diff_energies)
    plt_distri_rmse(path_plt, all_rmse)


def compare_scan_and_constrain_opt(result_path,tara_constrain_opt_path,xtb_sp_path):   # the main funtion
    # info = {}
    all_energies = {}
    #tera_infos = get_tera_mol(tara_scan_path)
    xtb_infos = get_xtb_energy_1(xtb_sp_path)
    #tera_infos = get_tera_from_opt(tara_opt_path)
    tera_infos = get_tera_constrain_energy(tara_constrain_opt_path)
    diff_energies = []
    all_rmse = []
    print(xtb_infos)
    print(tera_infos)
    for ik, xtb_value in xtb_infos.items():
        # print('start counting',ik,':')
        try:
            # diff_energies = []
            energies = []
            conf_rmse = []
            xtb_value_v = xtb_value.values()
            xtb_value_k = xtb_value.keys()
            length = len(xtb_value_k)
            tera_value = tera_infos[ik]

            for k,v in xtb_value.items():
                if k in tera_value.keys():
                    tera_value_v = tera_value[k]
                    en_xtb = v[0]
                    m_xtb = v[1]
                    en_tera = tera_value_v[0]
                    m_tera = tera_value_v[1]

                # print(Chem.MolToSmiles(m_xtb), Chem.MolToInchiKey(m_xtb))
                # print(Chem.MolToSmiles(m_tera), Chem.MolToInchiKey(m_tera))
                try:
                    rmse = AllChem.GetBestRMS(m_xtb, m_tera)
                except:
                    print('xtb:', Chem.MolToSmiles(m_xtb), Chem.MolToInchiKey(m_xtb))
                    print('tera:', Chem.MolToSmiles(m_tera), Chem.MolToInchiKey(m_tera))
                    print(ik, i, 'conformations change too much')
                    break

                diff_en = (float(en_xtb) - float(en_tera)) * 627.5094706
                diff_energies.append(diff_en)
                conf_rmse.append(rmse)
                all_rmse.append(rmse)
                energies.append([float(en_xtb), float(en_tera)])

                # diff_en = float(en_xtb) - float(en_tera)
                # diff_energies.append(diff_en)
                # conf_rmse.append(rmse)

                # energies.append([float(en_xtb), float(en_tera)])
            if len(energies) > 0:
                all_energies[ik] = energies
            path_plt = os.path.join(result_path, ik)
            if not os.path.exists((path_plt)):
                os.mkdir(path_plt)
            # plot_picture(path_plt, diff_energies, conf_rmse)
            plot_picture_1(path_plt, energies, conf_rmse)
            write_basic_info(path_plt, energies, conf_rmse)

        except:
            print(ik, 'is error')
            # info[ik] = [energies,conf_rmse]
    # '''
    # plt scatter
    path_plt = result_path
    # print(all_energies)
    plot_picture_2(path_plt, all_energies)
    plt_distri_energy(path_plt, diff_energies)
    plt_distri_rmse(path_plt, all_rmse)


if __name__ == "__main__":

    #analysis scan energies, rmse
    '''
    tara_scan_path = '/nfs2/zcc/new_comparing_process/frags/tera_scan'
    xtb_sp_path = '/nfs2/zcc/new_comparing_process/frags/xtb_scan_tera_sp'
    result_path = '/nfs2/zcc/new_comparing_process/frags/count_result_2'
    if not os.path.exists(result_path):
        os.mkdir(result_path)

    compare_rmsd_energy(result_path, tara_scan_path, xtb_sp_path)
    '''
     #analysis opt energies, rmse
    #'''
    result_path = '/nfs2/zcc/new_comparing_process/frags/count_result_vaccum_fix_opt'
    if not os.path.exists(result_path):
        os.mkdir(result_path)
    #tara_opt_path = '/nfs2/zcc/new_comparing_process/frag_xtb_opt_tera_strain'
    #xtb_sp_path ='/nfs2/zcc/new_comparing_process/frag_xtb_opt_tera_sp_finish/frag_xtb_opt_tera_sp'
    tera_opt_path = '/nfs2/zcc/new_comparing_process/frags/vaccum_xtb_opt_tera_strain_finish/vaccum_xtb_opt_tera_strain'
    xtb_sp_path = '/nfs2/zcc/new_comparing_process/frags/xtb_opt_tera_sp'
    compare_opt_energy_difference(result_path, tera_opt_path, xtb_sp_path)

    #'''
    '''
    # compare constrain scan and xtb_scan_tera_sp
    result_path = '/nfs2/zcc/new_comparing_process/count_result_constrain_scan'
    if not os.path.exists(result_path):
        os.mkdir(result_path)
    tara_constrain_opt_path = '/nfs2/zcc/new_comparing_process/frag_xtb_scan_tera_strain_finish/frag_xtb_scan_tera_strain'
    xtb_sp_path = '/nfs2/zcc/new_comparing_process/frag_xtb_scan_tera_sp_finish/frag_xtb_scan_tera_sp'
    compare_scan_and_constrain_opt(result_path,tara_constrain_opt_path,xtb_sp_path)
    '''






