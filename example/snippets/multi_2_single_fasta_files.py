import os
import sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
import utils
from utils import write2File
input_file= "/DATA/Dataset/BEPIPRED/BP3C50ID_external_test_set.fasta"
file_data = utils.multi_fasta_reader(input_file)
output_dir = "/DATA/EPITOPE/CLEANED_AACDB_DATA/dataset_features/BEPIPRED/fasta/"
for key, values in file_data.items():
    print(key)
    file_name = output_dir + str(key.replace(">","")) + ".fasta"
    utils.write2File(file_name, key+"\n"+values)
