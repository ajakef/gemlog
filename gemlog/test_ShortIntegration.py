from gemlog import *

def test_ReadGem_v0_9_single():
    ## This function accounts for most of the runtime, so it's a candidate for speeding up. Any flaws in its output (including hard-to-see flaws) can cause difficult-to-trace problems downstream.
    fn = 'data/FILE0001.077' 
    offset = 72000000.0 + 5263 # this is a reasonable offset to use; can't use 0 because bugs might slip through
    eps = 1e-12
    L1 = ReadGem_v0_9_single(fn, offset)
    L0 = slow_ReadGem_v0_9_single(fn, offset)
    assert np.abs(L1['data'] - L0['data']).max() < eps
    for i in ['gps', 'metadata']:
        for x in L0[i].keys():
            assert np.abs(L0[i][x] - L1[i][x]).max() < eps
            #assert L0[i]['t'].dtype == L1[i]['t'].dtype

def test_ReadGem():
    ## ReadGem always reads files in one block, so no sense in testing 25 files
    ReadGem(np.arange(3), 'data', SN = '077')


