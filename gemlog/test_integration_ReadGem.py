from gemlog import *

def test_ReadGem_empties():
    try:
        ReadGem(np.arange(1, 6), 'data', SN = '000') # test purely empty files
    except MissingRawFiles:
        print("Correctly catches missing raw files")
    except:
        raise Exception("Fails to catch missing raw files")
    else:
        raise Exception("Fails to catch missing raw files")

def test_ReadGem_edge_cases():
    try:
        ReadGem(np.array([23]), 'data', SN = '096') # test a malformed file
    except CorruptRawFile:
        print("Correctly catches corrupt raw file")
    except:
        raise Exception("Fails to catch corrupt raw file")
    else:
        raise Exception("Fails to catch corrupt raw file")

def test_ReadGem_good_data():
    ## ReadGem always reads files in one block, so no sense in testing 25 files
    ReadGem(np.arange(3), 'data', SN = '077') # test good data

