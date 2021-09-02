bash
Part A: run optmization
1. generate_xtb_tera.py:
to generate terachem and xtb optimization files
need open 'generate_xtb_tera.py' to modify input parameters
2. generate_tera_sh.py and generate_xtb_sh.py
to generate sh files
commond_format:python generate_tera_sh.py tera_opt/tera_scan<folder_name>
commond_format:python generate_xtb_sh.py xtb_opt/xtb_scan opt/scan
3. run_sbatch.py
commond_format:python run_sbatch.py tera_scan<folder_name> 0 30

Part B: run scan
4. generate_tera_scan.py generate_xtb_scan.py
to generate terachem and xtb scan files
need open scripts to modify input parameters
5. generate_tera_sh.py and generate_xtb_sh.py
6. run_sbatch.py

7.generate_xtb_sp_after_scanopt.py
# 'generate_start_xyz' method is use for terachem sp files based on xtb scan files
# 'generate_start_xyz_1' function is use for terachem sp files based on xtb opt files
need open the script and modify relative paths, then run it directly.
8. generate_tera_sh.py
commond format: python generate_tera_sh.py xtb_tera_sp<format_name>
9. run_sbatch.py

Part C: count result
10.count_cal_result.py
# 'compare_rmsd_energy' function is used for analyzing scan results
# 'compare_opt_energy_difference' is used for analyzing opt results

Part D:  generate_tera_freeze_opt.py
generate constrain optimizaion files based on xtb opt or scan
# 'generate_start_xyz_from_scan' function ..
# 'generate_start_xyz_from_opt' function ..

