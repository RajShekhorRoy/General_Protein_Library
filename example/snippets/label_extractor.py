import sys
import os
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
import utils
import os
# input_dir = "/DATA/AACDB/complex_structure/"
# interaction_files = "/DATA/AACDB/interacting_res_distance/"
# #Using heavy atom defination of 6A
# outputs_dir = "/DATA/EPITOPE/Output/"
# output_new_fasta_dir = outputs_dir+"fasta/"
# output_new_label_dir = outputs_dir+"labels/"
# output_new_structure_dir = outputs_dir+"clean_structure/"


def label_extractor_method(_input_struct_pdb,_interaction_values,_outputs_dir):
    # input_struct_dir= sys.argv[1]
    # interaction_values = sys.argv[2]
    # # Using heavy atom defination of 6A
    # outputs_dir = sys.argv[3]
    try: 
        input_struct_dir= _input_struct_pdb
        interaction_values=_interaction_values
        outputs_dir = _outputs_dir
        output_new_fasta_dir = outputs_dir + "fasta/"
        output_new_label_dir = outputs_dir + "labels/"
        antigen_structure_dir = outputs_dir + "antigen_structure/"
        antibody_structure_dir = outputs_dir + "antibody_structure/"
    
        #read interaction from their file
        # interact_data= utils.specific_filename_reader(interaction_files,".txt")
        # print(interact_data)
        
        
        
        print(interaction_values.split("_")[0],interaction_values.split("_")[1], interaction_values.split("_")[2])
        # find the pair of interaction
        # identify antigens
        combined_epitope_array=[]
        antigen_chain = interaction_values.split("_")[2]
        # identify antibody pair
        antibody_chains= interaction_values.split("_")[1]
        name = interaction_values.split("_")[0]
        complex_structure_file_name =input_struct_dir+name+"_"+antibody_chains+antigen_chain+".pdb"
        combnined_epitope_array=[]
        if os.path.exists(complex_structure_file_name):
            # print("Not found " + complex_structure_file_name)
            ###Remore Residue without CA
            complex_file  = utils.contents_to_info(utils.read_pdb(complex_structure_file_name))
            complex_file= utils.remove_residues_without_CA(complex_file)
            # print(complex_file)
        
        
            antigen_pdb = utils.separate_by_chain(complex_file, antigen_chain)
            antigen_fasta = utils.get_fasta_from_pdb_array(complex_file,antigen_chain)
            if len(antigen_fasta)>600:
               
                error_format = "{0}: {1}\n".format(_interaction_values, "great than max")
                utils.txt_appender(_filename="/DATA/EPITOPE/AACDB_data_processed/error.log", _text=error_format)
                raise Exception("great than max")
                return 
            fasta_file_name=output_new_fasta_dir+name+"_"+antigen_chain+".fasta"
            # print(antigen_fasta)
            fasta_str= utils.get_fasta_format(name+"_"+antigen_chain,antigen_fasta)
            utils.write2File(fasta_file_name, fasta_str)
            for antibody_chain in antibody_chains:
                # print(antigen_chain,antibody_chain)
                antibody_pdb = utils.separate_by_chain(complex_file, antibody_chain)
        
                antigen_serial_fix_a= utils.fix_res_num_atom(antigen_pdb)
                antibody_serial_fix_a = utils.fix_res_num_atom(antibody_pdb)
                utils.convert_to_pdb(antibody_serial_fix_a,  antibody_structure_dir+name+"_"+antibody_chain+".pdb")
                utils.convert_to_pdb(antigen_serial_fix_a, antigen_structure_dir+name+"_"+antigen_chain+".pdb")
                
                
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
                    # print((r + 1, c + 1, dist[r, c]))
                    combined_epitope_array.append(c + 1)
                utils.write2File(HRF_file_name, HRF_data_str)
            combined_labels_name = output_new_label_dir+interaction_values+".txt"
            combined_epitope_array= np.unique(combined_epitope_array)
            combined_epitope_str = utils.get_str_from_array(combined_epitope_array)
            utils.write2File(combined_labels_name, combined_epitope_str)
        
        # save combined epitope
    except Exception as e:
        print(e)
        error_format = "{0}: {1}\n".format(_interaction_values,e)
        utils.txt_appender(_filename="/DATA/EPITOPE/AACDB_data_processed/error.log", _text=error_format)
        
    
    return 
    
    
    


# clean the complexes
# save the new complex
#save the fastas
#save the labels