import utils

str_1 = ":A:B:C"
str_2 = ":A:B:C"

pdb_loc_1 = "./example/r-l_12_tidy.pdb"
pdb_loc_2 = "./example/r-l_23_tidy.pdb"
chain_map, overlapping_chain_1, overlapping_chain_2 = utils.chain_mapper(str_1, str_2)
pdb_1 = utils.contents_to_info(utils.read_pdb(pdb_loc_1))
pdb_2 = utils.contents_to_info(utils.read_pdb(pdb_loc_2))
# print(pdb_1)
#### use mapped chains

fasta_1 = utils.get_fasta_from_pdb_array(pdb_1)
fasta_2 = utils.get_fasta_from_pdb_array(pdb_2)

chains_1 = utils.get_unique_chains(pdb_1)
chains_2 = utils.get_unique_chains(pdb_2)
iq_profile_obj = utils.iq_profile()

##Inital data collection part

for chain in chains_1:
    an_original_pdb = utils.separate_by_chain(pdb_1, chain)
    a_fasta = utils.get_fasta_from_pdb_array(pdb_1)
    iq_profile_obj.original_fasta_1[chain] = a_fasta
    iq_profile_obj.original_pdb_1[chain] = utils.fix_serial(an_original_pdb)
    iq_profile_obj.overlapping_pdb_1_ca[chain] = utils.get_skeleton(an_original_pdb)

for chain in chains_2:
    an_original_pdb = utils.separate_by_chain(pdb_2, chain)
    a_fasta = utils.get_fasta_from_pdb_array(pdb_2)
    iq_profile_obj.original_fasta_2[chain] = a_fasta
    iq_profile_obj.original_pdb_2[chain] = utils.fix_serial(an_original_pdb)
    iq_profile_obj.overlapping_pdb_2_ca[chain] = utils.get_skeleton(an_original_pdb)

iq_profile_obj.original_interactions = utils.get_all_interactions(chains_1, iq_profile_obj)
print(iq_profile_obj.original_interactions)

#get overlapping
for a_ref_chains in iq_profile_obj.original_interactions:
    print(a_ref_chains)
    print(chain_map.get(a_ref_chains))
    fasta_1 = iq_profile_obj.original_fasta_1[a_ref_chains]
    fasta_2 = iq_profile_obj.original_fasta_2[chain_map.get(a_ref_chains)]
    utils.find_common_fasta(fasta_1,fasta_2)

#calculate i_rmsd after cleaning

