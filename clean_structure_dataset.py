import numpy as np

import utils
import os
input_dir = "/Users/rajshekhorroy/Documents/small_data/complex_structure/"
interaction_files = "/Users/rajshekhorroy/Documents/small_data/interacting_res_distance/"
#Using heavy atom defination of 6A
outputs_dir = "/Users/rajshekhorroy/Documents/output/"
output_new_fasta_dir = outputs_dir+"fasta/"
output_new_label_dir = outputs_dir+"labels/"
output_new_structure_dir = outputs_dir+"clean_structure/"


#read interaction from their file
interact_data= utils.specific_filename_reader(interaction_files,".txt")
# print(interact_data)
for values in interact_data:
    print(values.split("_")[0],values.split("_")[1], values.split("_")[2])
    # find the pair of interaction
    # identify antigens
    combined_epitope_array=[]
    antigen_chain = values.split("_")[2]
    # identify antibody pair
    antibody_chains= values.split("_")[1]
    name = values.split("_")[0]
    complex_structure_file_name =input_dir+name+"_"+antibody_chains+antigen_chain+".pdb"
    combnined_epitope_array=[]
    if os.path.exists(complex_structure_file_name):
        # print("Not found " + complex_structure_file_name)
        ###Remore Residue without CA
        complex_file  = utils.contents_to_info(utils.read_pdb(complex_structure_file_name))
        complex_file= utils.remove_residues_without_CA(complex_file)
        # print(complex_file)


        antigen_pdb = utils.separate_by_chain(complex_file, antigen_chain)
        antigen_fasta = utils.get_fasta_from_pdb_array(complex_file,antigen_chain)
        fasta_file_name=output_new_fasta_dir+name+"_"+antigen_chain+".fasta"
        print(antigen_fasta)
        fasta_str= utils.get_fasta_format(name+"_"+antigen_chain,antigen_fasta)
        utils.write2File(fasta_file_name, fasta_str)
        for antibody_chain in antibody_chains:
            print(antigen_chain,antibody_chain)
            antibody_pdb = utils.separate_by_chain(complex_file, antibody_chain)

            antigen_serial_fix_a= utils.fix_res_num_atom(antigen_pdb)
            antibody_serial_fix_a = utils.fix_res_num_atom(antibody_pdb)
            if name =="1DZB":
                print("wait here")

            try:

                dist,contact = utils.get_distance_map_heavy_atom(antibody_serial_fix_a,antigen_serial_fix_a)
                dist_file_name = output_new_label_dir+name+"_dist_"+antigen_chain+antibody_chain+".txt"
                contact_file_name = output_new_label_dir + name + "_contact_6_" + antigen_chain + antibody_chain + ".txt"
                HRF_file_name = output_new_label_dir + name + "_HRF_6_" + antigen_chain + antibody_chain + ".txt"

                np.savetxt(dist_file_name, dist, fmt="%.3f")
                np.savetxt(contact_file_name, contact, fmt="%.3f")
                ## contacts and distance in term of interace
                HRF_data_str=""
                nz = np.argwhere(dist <= 6)
                #save contacts in HRF
            except Exception as  e:
                print( "something went wrong "+ e )


            for r, c in nz:
                HRF_data_str += "{0},{1},{2}".format(r + 1, c + 1, dist[r, c])+"\n"
                print((r + 1, c + 1, dist[r, c]))
                combined_epitope_array.append(c + 1)
            utils.write2File(HRF_file_name, HRF_data_str)
        combined_labels_name = output_new_label_dir+values+".txt"
        combined_epitope_array= np.unique(combined_epitope_array)
        combined_epitope_str = utils.get_str_from_array(combined_epitope_array)
        utils.write2File(combined_labels_name, combined_epitope_str)

        # save combined epitope







# clean the complexes
# save the new complex
#save the fastas
#save the labels