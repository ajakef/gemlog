from gemlog import *

def test_ReadGem_missing():
    try:
        ReadGem(np.arange(10, 20), 'data', SN = '000') # test missing files
    except MissingRawFiles:
        print("Correctly catches missing raw files")
    except:
        raise Exception("Fails to catch missing raw files")
    else:
        raise Exception("Fails to catch missing raw files")

def test_ReadGem_missing():
    try:
        ReadGem(np.arange(5), 'data', SN = '000') # test purely empty files
    except MissingRawFiles:
        print("Correctly catches empty/missing raw files")
    except:
        raise Exception("Fails to catch empty/missing raw files")
    else:
        raise Exception("Fails to catch empty/missing raw files")

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


## Convert tests: ensure that it doesn't crash, and that the output mseed file is identical to a reference
def test_Convert_good_data():
    Convert('data', SN = '077', convertedpath = 'test_output_mseed')
    output = obspy.read('test_output_mseed/2020-04-24T22:00:00..077..HDF.mseed')
    reference = obspy.read('data/2020-04-24T22:00:00..077..HDF.mseed')
    assert output.__eq__(reference)
