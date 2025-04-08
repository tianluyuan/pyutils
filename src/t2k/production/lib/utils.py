import pickle

# this overrides the baseline geometery setting and can be
# used for producing neutrino interactions on a specific target
STANDALONE_GEOMETRIES = ('TPC', 'FGD', 'P0D', 'DSECAL')


def parseMissingFileList(aFile, specifier='missing'):
    with open(aFile) as f:
        missingDict = pickle.load(f)

    return missingDict[specifier]


def get_run_and_subrun(filename):
    dash = filename.find('-')
    run = filename[dash-8:dash]
    subrun = filename[dash+1:dash+5]
    return int(run), int(subrun)
