from gemlog import *

def test_ReadGem_empties():
    ReadGem(np.arange(1, 6), 'data', SN = '000') # test purely empty files

def test_ReadGem_edge_cases():
    ReadGem(np.array([23]), 'data', SN = '096') # test a malformed file

def test_ReadGem_good_data():
    ## ReadGem always reads files in one block, so no sense in testing 25 files
    ReadGem(np.arange(3), 'data', SN = '077') # test good data

