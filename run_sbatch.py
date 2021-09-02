import os
import sys


def run_qm_scripts(input_folder, init_p,end_p):
    init_p = int(init_p)
    end_p = int(end_p)
    chmod_command = "chmod -R 777 %s"%input_folder
    os.system(chmod_command) 
    sub_dir_list = os.listdir(input_folder)
    sec_para = len(sub_dir_list)
    if end_p <= len(sub_dir_list):
        sec_para = end_p
    
    for inchi_key in sub_dir_list[init_p:sec_para]:
        full_dir_path = os.path.join(input_folder, inchi_key)
        all_file_list = os.listdir(full_dir_path)
        os.chdir(full_dir_path)
        #print(os.chdir(full_dir_path))
        if len(all_file_list) > 0:
            for file in all_file_list:
                if file.endswith(".sh"):
                    file_path = os.path.join(full_dir_path,file)
                    #print("sbatch %s"%file_path)
                    os.system("sbatch %s"%file_path)

    record_txt = input_folder+"_record.txt"
    with open(record_txt,"a")as f:
        f.write("write from {} to {}\n".format(init_p,sec_para))
        f.write("all the number of jobs are {}\n".format(len(sub_dir_list)))
        for inchi_key in sub_dir_list[init_p:sec_para]:
            f.write(inchi_key+"\n")
        f.close()

if __name__ == "__main__":
    curr_path = os.getcwd()
    input_folder = os.path.join(curr_path, sys.argv[1])
    init_p = sys.argv[2]
    end_p = sys.argv[3]
    #output_folder = input_folder+"_out"
    run_qm_scripts(input_folder=input_folder,init_p=init_p, end_p=end_p)

