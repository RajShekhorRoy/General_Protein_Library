
import utils
str_1=":A:B:C"
str_2=":A:B:C"

pdb_loc_1= "./example/r-l_12_tidy.pdb"
pdb_loc_2= "./example/r-l_23_tidy.pdb"

pdb_1= utils.contents_to_info(utils.read_pdb(pdb_loc_1))
pdb_2= utils.contents_to_info(utils.read_pdb(pdb_loc_2))
print(pdb_1)
#### use mapped chains

fasta_1=utils.get_fasta_from_pdb_array(pdb_1)
fasta_2=utils.get_fasta_from_pdb_array(pdb_2)

chains_1 = utils.get_unique_chains(pdb_1)
chains_2 = utils.get_unique_chains(pdb_2)
iq_profile_obj=  utils.iq_profile()
for chain in chains_1:
    an_original_pdb = utils.separate_by_chain(pdb_1, chain)
    a_fasta = utils.get_fasta_from_pdb_array(pdb_1)

    iq_profile_obj.original_fasta_1[chain]=a_fasta

    iq_profile_obj.original_pdb_1[chain]=an_original_pdb



for chain in chains_2:
    an_original_pdb = utils.separate_by_chain(pdb_2, chain)
    a_fasta = utils.get_fasta_from_pdb_array(pdb_2)

    iq_profile_obj.original_fasta_2[chain]=a_fasta

    iq_profile_obj.original_pdb_2[chain]=an_original_pdb

print(iq_profile_obj)