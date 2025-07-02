# pdb_data = list(filter(lambda x: ( x.alt_loc != 'B'), copy.deepcopy(pdb_data)))

import utils
pdb_loc_1 = "../r-l_12_tidy_unclean.pdb"
out_loc = "../r-l_12_tidy_clean_A.pdb"
pdb_1 = utils.contents_to_info(utils.read_pdb(pdb_loc_1))
# print(pdb_1)
chains = utils.get_unique_chains(pdb_1)
# pdb_data = list(filter(lambda x: ( x.alt_loc != 'B'), copy.deepcopy(pdb_data)))
import copy
import utils
pdb_loc_1 = "../r-l_12_tidy_unclean.pdb"
out_loc = "../r-l_12_tidy_clean_All.pdb"
pdb_1 = utils.contents_to_info(utils.read_pdb(pdb_loc_1))
# print(pdb_1)
chains = utils.get_unique_chains(pdb_1)
all_pdb = []
for values in chains:
    fix_res_num_pdb =utils.fix_res_num_atom(utils.separate_by_chain(pdb_1,values))
    alt_loc_removed_pdb = list(filter(lambda x: ( x.alt_loc != 'B'), copy.deepcopy(fix_res_num_pdb)))
    # atom_num_fixed = utils.fix_serial(alt_loc_removed_pdb)
    all_pdb.extend(alt_loc_removed_pdb)

cleaned_pdb= utils.pdb_from_array_multi_chain(all_pdb,out_loc)