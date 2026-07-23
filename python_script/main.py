import subprocess
import sys
import os
import pandas as pd 
from contextlib import redirect_stdout, redirect_stderr
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed

## DEFINITION OF SCRIPT PATH ##
theo_SLiM = "/mnt/DATA1/SCRIPT_SLIM/theoritical_simulations.slim"
life_cycle_Dsarg_SLiM = "/mnt/DATA1/SCRIPT_SLIM/life_cycle_Dsarg_simulations.slim"
life_cycle_Emarg_SLiM = "/mnt/DATA1/SCRIPT_SLIM/life_cycle_Emarg_simulations.slim"
laplace_law_SLiM = "/mnt/DATA1/SCRIPT_SLIM/laplace_law_simulations.slim"
student_law_SLiM = "/mnt/DATA1/SCRIPT_SLIM/student_law_simulations.slim"
herm_Dsarg_SLiM_500 = "/mnt/DATA1/SCRIPT_SLIM/hermaphrodism_Dsarg_500_simulations.slim"
herm_Dsarg_SLiM_400 = "/mnt/DATA1/SCRIPT_SLIM/hermaphrodism_Dsarg_400_simulations.slim"
herm_Dsarg_SLiM_300 = "/mnt/DATA1/SCRIPT_SLIM/hermaphrodism_Dsarg_300_simulations.slim"
herm_Emarg_SLiM_5000 = "/mnt/DATA1/SCRIPT_SLIM/hermaphrodism_Emarg_5000_simulations.slim"
herm_Emarg_SLiM_4000 = "/mnt/DATA1/SCRIPT_SLIM/hermaphrodism_Emarg_4000_simulations.slim"
herm_Emarg_SLiM_3000 = "/mnt/DATA1/SCRIPT_SLIM/hermaphrodism_Emarg_3000_simulations.slim"
hetero_m_SLiM = "/mnt/DATA1/SCRIPT_SLIM/heterogeneous_m_simulations_v2.slim"
hetero_N_SLiM = "/mnt/DATA1/SCRIPT_SLIM/heterogeneous_N_simulations_v2.slim"
grad_m_SLiM = "/mnt/DATA1/SCRIPT_SLIM/gradual_m_simulations_v2.slim"
grad_N_SLiM = "/mnt/DATA1/SCRIPT_SLIM/gradual_N_simulations_v2.slim"
temporal_m_SLiM = "/mnt/DATA1/SCRIPT_SLIM/temporal_m_simulations.slim"
temporal_N_SLiM = "/mnt/DATA1/SCRIPT_SLIM/temporal_N_simulations.slim"
barrier_m_SLiM = "/mnt/DATA1/SCRIPT_SLIM/barrier_m_simulations.slim"
barrier_m_v2_SLiM = "/mnt/DATA1/SCRIPT_SLIM/barrier_m_simulations_v2.slim"

### DEFINITION OF FUNCTION ###
# Function to execute a Python script
def run_python_script(script_path, args=None, conda_env=None):
   try: 
       if conda_env: 
           command = ["conda", "run", "-n", conda_env, "python", script_path]
       else:
           command = [sys.executable, script_path]  
       if args:
           command.extend(map(str, args))
       print(f"\n>>> Running {script_path} with args {args} in env {conda_env}")
       subprocess.run(command, check=True, stderr=subprocess.STDOUT)
   except subprocess.CalledProcessError as e:
       print(f"❌ Error executing {script_path} : {e}")
       sys.exit(1)

# Function to execute a SLiM script
def run_slim_script(slim_script_path, args=None):
   try:
       command = ["slim"]
       if args:
           command.extend(map(str, args))
       command.append(slim_script_path) 
       print(f"\n>>> Running SLiM script {slim_script_path} with args {args}")
       subprocess.run(command, check=True, stderr=subprocess.STDOUT)
   except subprocess.CalledProcessError as e:
       print(f"❌ Error executing SLiM script : {e}")
       sys.exit(1)

# Function to run a single simulation 
def run_single_simulation(row_dict):
    # Force all intern call to use only 1 threads
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"

    # Create a log file for each simulation
    log_file_name = f"logs/sim_{row_dict['name']}.log" 
    os.makedirs("logs", exist_ok=True)

    with open(log_file_name, "w") as log_file, redirect_stdout(log_file), redirect_stderr(log_file):
        try:
            print(f"\n{'='*30}\n▶ Step 1: Parsing parameters\n{'='*30}")
            print(f"🔧 Starting simulation for: {row_dict['name']}")
       
            ## Extract parameters ## 
            name = row_dict['name']
            spatial = row_dict['spatial']
            m = row_dict['m']
            N = row_dict['N']
            dispersal_law = row_dict['dispersal_law']
            life_cycle = row_dict['life_cycle']
            herm = row_dict['hermaphrodism']
            K = row_dict['K']
            t = row_dict['t']
            rep = row_dict['rep']
            model = row_dict['model']
            
            if model == "life_cycle_Dsarg":
                slim_script_to_run = life_cycle_Dsarg_SLiM
            elif model == "life_cycle_Emarg":
                slim_script_to_run = life_cycle_Emarg_SLiM
            elif model == "herm_Dsarg_500":
                slim_script_to_run = herm_Dsarg_SLiM_500
            elif model == "herm_Dsarg_400":
                slim_script_to_run = herm_Dsarg_SLiM_400
            elif model == "herm_Dsarg_300":
                slim_script_to_run = herm_Dsarg_SLiM_300
            elif model == "herm_Emarg_5000":
                slim_script_to_run = herm_Emarg_SLiM_5000
            elif model == "herm_Emarg_4000":
                slim_script_to_run = herm_Emarg_SLiM_4000
            elif model == "herm_Emarg_3000":
                slim_script_to_run = herm_Emarg_SLiM_3000    
            elif dispersal_law == "Laplace":
                slim_script_to_run = laplace_law_SLiM
            elif dispersal_law == "Student":
                slim_script_to_run = student_law_SLiM
            elif spatial == "hetero_N":
                slim_script_to_run = hetero_N_SLiM
            elif spatial == "hetero_m":
                slim_script_to_run = hetero_m_SLiM
            elif spatial == "grad_N":
                slim_script_to_run = grad_N_SLiM
            elif spatial == "grad_m":
                slim_script_to_run = grad_m_SLiM
            elif spatial == "barrier_m":
                slim_script_to_run = barrier_m_SLiM
            elif model == "temporal_m":
                slim_script_to_run = temporal_m_SLiM
            elif model == "temporal_N":
                slim_script_to_run = temporal_N_SLiM
            elif spatial == "barrier_m_v2":
                slim_script_to_run = barrier_m_v2_SLiM
            else:
                slim_script_to_run = theo_SLiM
            
            ## Define output folder ##
            output_TS = f"/mnt/DATA1/TS/"
            output_IBD = f"/mnt/DATA1/IBD/"
            output_IBDAnalysis = f"/mnt/DATA1/INFERENCE/IBD_ANALYSIS/"
            output_MAPS = f"/mnt/DATA1/INFERENCE/MAPS/"
            map_file = f"/mnt/DATA1/MAP/map_100x100_grayscale.png"
            
            ## Obtain IBD segments ## 
            # Step 0 : Run SLiM script 
            print(f"\n{'='*30}\n▶ Step 2: Running SLiM\n{'='*30}")
            ts_file = f"{output_TS}forw_{spatial}_m{m}_N{N}_{model}_rep{rep}.trees"
            if t == 201 and K == 200:
                if not os.path.exists(ts_file):
                    print(f"Running SLiM because output file not found: {ts_file}")
                    run_slim_script(slim_script_to_run, args=["-d", f"S_disp={m}", "-d", f"K={N}", "-d", f"rep={rep}", "-d", f'model="{model}"', "-d", f'spatial="{spatial}"'])
                else:
                    print(f"Skipping SLiM (file already exists): {ts_file}")
            else:
                if t != 201:
                    print(f"Skipping SLiM because t = {t} != 201 (using theo TS file)")
                if K != 200:
                    print(f"Skipping SLiM because K = {K} != 200 (using theo TS file)")
                if not os.path.exists(ts_file):
                    raise FileNotFoundError(f"Expected TS file not found: {ts_file}")

            # Step 1 : Extract IBD from tree sequence file 
            print(f"\n{'='*30}\n▶ Step 3: Extracting IBD segment\n{'='*30}")
            ibd_file = f"{output_IBD}IBD_coord_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            if not os.path.exists(ibd_file):
                run_python_script(TS2IBD, args=["-i", ts_file, "-o", output_IBD, "-hap", str(K), "-t_max", str(t), '-loc'], conda_env = "ts-env")
            else:
                print(f"Skipping IBD Extraction (file already exists): {ibd_file}")

            ## Plot ## 
            # Step 2 : Plot the length distribution of IBD segments 
            print(f"\n{'='*30}\n▶ Step 4: Plotting IBD segments statistics\n{'='*30}")
            ibd_file = f"{output_IBD}IBD_coord_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            distrib_file = f"{output_IBD}distrib_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.png"
            if not os.path.exists(distrib_file):
                run_python_script(distrib, args=["-i", ibd_file, "-o", output_IBD, "-c", "-hap", str(K)])
            else:
                print(f"Skipping Distribution Plot (file already exists): {distrib_file}")

            # Step 3 : Plot the genetic distance by geographic distance 
            loc_file = f"{output_IBD}Ind_Coord_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            ibdist_file = f"{output_IBD}Prop_vs_distGeo_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.png"
            if not os.path.exists(ibdist_file):
                run_python_script(isolation_by_distance, args=["-i", ibd_file, "-o", output_IBD, "-loc", loc_file, "-d"], conda_env = "ts-env")
            else:
                print(f"Skipping IBD-Distance (file already exists): {ibdist_file}")

            # Step 4 : Plot the sampled individuals on the map 
            coord_file = f"{output_IBD}Coord_all_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            coord_sampled_file = f"{output_IBD}Coord_samp_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            map_indivs_file = f"{output_IBD}Map_individuals_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.png"
            if not os.path.exists(map_indivs_file):
                run_python_script(map_location, args=["-i", ibd_file, "-o", output_IBD, "-coord_sim", coord_file, "-coord", coord_sampled_file, "-map", map_file], conda_env = "ts-env")
            else:
                print(f"Skipping Map Plot (file already exists): {map_indivs_file}")

            ## IBD-Analysis ## 
            print(f"\n{'='*30}\n▶ Step 5: Running IBD-Analysis\n{'='*30}")
            # Step 5 : Convert IBD file into IBD-Analysis file 
            IBDAnalysis_ibd_file = f"{output_IBD}IBD-Analysis_IBD_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            if not os.path.exists(IBDAnalysis_ibd_file):
                run_python_script(TOIBD_Analysis, args=["-i", ibd_file, "-o", output_IBD, "-loc", loc_file])
            else:
                print(f"Skipping IBD-Analysis Conversion (file already exists): {IBDAnalysis_ibd_file}")

            # Step 6 : Run IBD-Analysis 
            IBDAnalysis_ibd_file = f"{output_IBD}IBD-Analysis_IBD_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            IBDAnalysis_coord_file = f"{output_IBD}IBD-Analysis_Coord_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            IBDAnalysis_results = f"{output_IBDAnalysis}Results_MLE_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            if not os.path.exists(IBDAnalysis_results):
                run_python_script(MLE_IBD_Analysis, args=["-i", IBDAnalysis_ibd_file, "-o", output_IBDAnalysis, "-loc", IBDAnalysis_coord_file, "-p", str(m), str(N)]) 
            else:
                print(f"Skipping IBD-Analysis Run (file already exists): {IBDAnalysis_results}")

            ## MAPS ##
            print(f"\n{'='*30}\n▶ Step 6: Running MAPS\n{'='*30}")
            # Step 7 : Convert IBD and coordinates file into MAPS file 
            output_MAPS_file = f"{output_MAPS}{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}/"
            os.makedirs(output_MAPS_file, exist_ok=True)
            MAPS_file = f"{output_MAPS_file}{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}_2_Inf.sims"
            if not os.path.exists(MAPS_file):
                run_python_script(TOMAPS, args=["-i", ibd_file, "-o", output_MAPS_file, "-loc", loc_file])
            else: 
                print(f"Skipping MAPS Conversion (file already exists): {MAPS_file}")
            
            # Step 8 : Create parameters file  
            data_path = f"{output_MAPS_file}{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}_2_Inf"
            mcmc_path = f"{output_MAPS_file}2_Inf"
            params_name = f"params_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.ini"
            params_file = f"{output_MAPS}params_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.ini"
            haplotypes = K
            if not os.path.exists(params_file):
                run_python_script(params, args=["-o", output_MAPS, "-f", params_name, "-i", data_path, "-m", mcmc_path, "-hap", haplotypes])
            else:
                print(f"Skipping Parameters File Creation (file already exists): {params_file}")
            
            # Step 9 : Run MAPS 
            results_folder = f"{output_MAPS_file}2_Inf"
            os.makedirs(results_folder, exist_ok=True)
            params_file = f"{output_MAPS}params_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.ini"
            results_file = f"{output_MAPS}Results_MAPS_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.txt"
            if not os.path.exists(results_file):
                run_python_script(MAPS, args=["-i", output_MAPS_file, "-p", params_file])
            else:
                print(f"Skipping MAPS Run (file already exists): {results_file}")
            
            # Step 10 : Plot the results   
            summary_file = f"{output_MAPS}Summary_{spatial}_m{m}_N{N}_{model}_K{K}_t{t}_rep{rep}.csv"
            if not os.path.exists(summary_file):
                run_python_script(plotmaps, args=["-i", results_folder, "-o", output_MAPS])
            else:
                print(f"Skipping Results Plot (file already exists): {summary_file}") 

            print(f"\n{'='*30}\n▶ Step 6: Finished Simulation\n{'='*30}")
            return f"✅ SUCCESS : {row_dict['name']}"
        
        except Exception as e:
            error_msg = f"\n❌ ERROR in {row_dict['name']} : {str(e)}\n"
            print(error_msg)
            with open(log_file_name, "a") as log_file:
                log_file.write(error_msg)
            return f"❌ ERROR : {row_dict['name']} : {str(e)}"

    
## RUN IN PARALLEL ##
if __name__ == "__main__":
    df = pd.read_csv("/mnt/DATA1/Parameters.csv")

    # Convert rows to dict to pass to subprocesses (safer que passer Series)
    row_dicts = [row.to_dict() for _, row in df.iterrows()]
    
    # Run 30 simulations at same time 
    max_parallel_runs = 30 
    results = []
    
    with ProcessPoolExecutor(max_workers=max_parallel_runs) as executor:
        futures = {executor.submit(run_single_simulation, row): row for row in row_dicts}
        for future in as_completed(futures):
            try:
                result = future.result()
                print(result)  # Print SUCCESS or ERROR message
                results.append(result)
            except Exception as e:
                print(f"❌ Unexpected crash in one of the simulations: {str(e)}")
                
    # All simulations finished 
    print("\n>>> All simulations finished : Concatenate results\n")

    MAPS_folder = "/mnt/DATA1/INFERENCE/MAPS"
    IBDAnalysis_folder = "/mnt/DATA1/INFERENCE/IBD_ANALYSIS"
    output_folder = "/mnt/DATA1"
    run_python_script(concat, args=["-M", MAPS_folder, "-I", IBDAnalysis_folder, "-o", output_folder])

