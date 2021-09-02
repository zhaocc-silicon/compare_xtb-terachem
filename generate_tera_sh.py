import os
import sys


def create_qm_scripts(input_folder,  node=""):

    l0 = "#!/bin/bash\n"
    #	command = RunCommand()
    sub_dir_list = os.listdir(input_folder)
    for inchi_key in sub_dir_list:
        full_dir_path = os.path.join(input_folder, inchi_key)
        all_file_list = os.listdir(full_dir_path)

        t_file_list = []
        if len(all_file_list)>0:
            for file in all_file_list:
                if file.endswith(".start"):
                    t_file_list.append(file)

        if len(t_file_list)>0:
            for t_file in t_file_list:
                t_name = t_file.replace(".start", "")
                out_name = t_file.replace(".start", ".out")
                full_sh_path = os.path.join(full_dir_path, "".join([t_name, ".sh"]))
                #gjf_path = os.path.join(full_dir_path, t_file)

                l1 = "#SBATCH --nodes=1\n"
                l2 = "#SBATCH --oversubscribe\n"
                l2_2 = "#SBATCH --gres=gpu:1\n"
                l3 = "#SBATCH --ntasks=1\n"
                l4 = "#SBATCH --cpus-per-task=2\n"
                l5 = "#SBATCH --partition=project\n"
                l6 = "#SBATCH --qos=maxjobs\n"
                #l7 = "#SBATCH --time=10:00:00\n"
                l8 = "#SBATCH --nice=1\n"
                l9 = "source /software/TeraChem_1.94-beta/SetTCVars.sh\n"
                l10 = "terachem  {}   >{}\n".format(t_file,out_name)

                with open(full_sh_path, "w")as f:  # open sh file
                    f.write(l0)
                    f.write(l1)
                    f.write(l2)
                    f.write(l2_2)
                    f.write(l3)
                    f.write(l4)
                    f.write(l5)
                    f.write(l6)
                    #f.write(l7)
                    f.write(l8)
                    f.write(l9)
                    f.write(l10)
                    f.close()
if __name__ == "__main__":
    curr_path = os.getcwd()
    input_folder = os.path.join(curr_path, sys.argv[1])
    #output_folder = input_folder+"_out"
    create_qm_scripts(input_folder=input_folder)
