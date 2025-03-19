import utils

str_1 = ":A:B:C"
str_2 = ":A:B:C"

pdb_loc_1 = "./example/r-l_12_tidy.pdb"
pdb_loc_2 = "./example/r-l_23_tidy.pdb"
##read and clean Pdb
chain_map, overlapping_chain_1, overlapping_chain_2 = utils.chain_mapper(str_1, str_2)
pdb_1 = utils.contents_to_info(utils.read_pdb(pdb_loc_1))
pdb_2 = utils.contents_to_info(utils.read_pdb(pdb_loc_2))
# print(pdb_1)
##maybe common chain mappings
#introduce mask if missing


#### use mapped chains

fasta_ref = utils.get_fasta_from_pdb_array(pdb_1)
fasta_com = utils.get_fasta_from_pdb_array(pdb_2)

chains_ref = utils.get_unique_chains(pdb_1)
chains_com = utils.get_unique_chains(pdb_2)
iq_profile_obj = utils.iq_profile()

##Inital data collection part

for chain in chains_ref:
    an_original_pdb = utils.separate_by_chain(pdb_1, chain)
    a_fasta = utils.get_fasta_from_pdb_array(pdb_1)
    iq_profile_obj.original_ref_fasta[chain] = a_fasta
    iq_profile_obj.original_ref_pdb[chain] = utils.fix_serial(an_original_pdb)
    iq_profile_obj.original_ref_pdb_ca[chain] = utils.get_skeleton(an_original_pdb)

for chain in chains_com:
    an_original_pdb = utils.separate_by_chain(pdb_2, chain)
    a_fasta = utils.get_fasta_from_pdb_array(pdb_2)
    iq_profile_obj.original_com_fasta[chain] = a_fasta
    iq_profile_obj.original_com_pdb[chain] = utils.fix_serial(an_original_pdb)
    iq_profile_obj.original_com_pdb_ca[chain] = utils.get_skeleton(an_original_pdb)

iq_profile_obj.ref_interactions = utils.get_all_interactions(chains_ref, iq_profile_obj)
iq_profile_obj.com_interactions = utils.get_all_interactions(chains_com, iq_profile_obj)
# print(iq_profile_obj.original_interactions)

#get overlapping
for a_ref_chains in iq_profile_obj.ref_interactions:
    print(a_ref_chains)
    print(chain_map.get(a_ref_chains))
    fasta_1 = iq_profile_obj.original_ref_fasta[a_ref_chains]
    fasta_2 = iq_profile_obj.original_com_fasta[chain_map.get(a_ref_chains)]
    iq_profile_obj.aligned_ref_fasta[a_ref_chains], iq_profile_obj.aligned_com_fasta[a_ref_chains], iq_profile_obj.overlapping_fastas[a_ref_chains] =  utils.find_common_fasta(fasta_1,fasta_2)
# calculate i_rmsd after cleaning
#get contact map
for a_ref_chains in iq_profile_obj.ref_interactions:
    for chain in iq_profile_obj.ref_interactions[a_ref_chains]:
        print(a_ref_chains,chain)
        dist_map,contact_map = utils.get_distance_map(first_chain_CA=iq_profile_obj.original_ref_pdb_ca[a_ref_chains],second_chain_CA=iq_profile_obj.original_ref_pdb_ca[chain])
        iq_profile_obj.distance_maps_ref[a_ref_chains+"_"+chain],     iq_profile_obj.contact_maps_ref[a_ref_chains+"_"+chain] = dist_map,contact_map


for a_com_chains in iq_profile_obj.com_interactions:
    for chain in iq_profile_obj.com_interactions[a_com_chains]:
        print(a_com_chains,chain)
        dist_map ,contact_map= utils.get_distance_map(first_chain_CA=iq_profile_obj.original_com_pdb_ca[a_com_chains],second_chain_CA=iq_profile_obj.original_com_pdb_ca[chain])
        iq_profile_obj.distance_maps_com[a_com_chains+"_"+chain] ,     iq_profile_obj.contact_maps_com[a_com_chains+"_"+chain]= dist_map,contact_map
