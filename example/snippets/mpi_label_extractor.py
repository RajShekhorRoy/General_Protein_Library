import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
import utils
from label_extractor import label_extractor_method


interaction_files = "/DATA/AACDB/interacting_res_distance_/"
interact_data= utils.specific_filename_reader(interaction_files,".txt")
input_dir = "/DATA/AACDB/complex_structure/"
outputs_dir = "___/DATA/EPITOPE/removed_alt_loc___/"

 
#_input_struct_pdb,_interaction_values,_outputs_dir
job_array= []
for values in interact_data:
    job_array.append((input_dir,values,outputs_dir))


max_parallel =15

with ProcessPoolExecutor(max_workers=max_parallel) as executor:
    futures = [executor.submit(label_extractor_method , *t ) for t in job_array]

    for future in as_completed(futures):
        print(future.result())  # printed as so