import utils
import os
input_dir = "/DATA/AACDB/complex_structure/"
outputs_dir = "/DATA/EPITOPE/test_output/pdb/"
output_new_fasta_dir = "/DATA/EPITOPE/test_output/fasta/"
output_new_label_dir = "/DATA/EPITOPE/test_output/labels/"

interaction_files = "/DATA/AACDB/interacting_res_distance/"

#read interaction from their file
interact_data= utils.specific_filename_reader(interaction_files,".txt")[0:10]
# print(interact_data)
for values in interact_data:
    print(values.split("_")[0],values.split("_")[1], values.split("_")[2])
    # find the pair of interaction
    # identify antigens
    antigen_chain = values.split("_")[2]
    # identify antibody pair
    antibody_chain= values.split("_")[1]
    name = values.split("_")[0]
    complex_structure_file_name =input_dir+name+"_"+antibody_chain+antigen_chain+".pdb"
    
    if os.path.exists(complex_structure_file_name):
        print("Not found " + complex_structure_file_name)

        complex_file  = utils.read_pdb(complex_structure_file_name)
        print(complex_file)
 
        an_original_pdb = utils.separate_by_chain(pdb_2, chain)
        a_fasta = utils.get_fasta_from_pdb_array(pdb_2, chain)
        iq_profile_obj.original_com_fasta[chain] = a_fasta
        iq_profile_obj.original_com_pdb[chain] = utils.fix_serial(an_original_pdb)
        iq_profile_obj.original_com_pdb_ca[chain] = utils.get_skeleton(an_original_pdb)
        output_name = outputs_dir+name+"_"+antibody_chain+antigen_chain+".pdb"
    print(f"Saved renumbered PDB as {output_name}")




# clean the complexes
# save the new complex
#save the fastas
#save the labels