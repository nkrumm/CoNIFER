import statsmodels.api as sm
import argparse
import cPickle
import tables
import pandas
import glob
import os

data = sm.datasets.scotland.load()
data.exog = sm.add_constant(data.exog)

def loadRPKM(rpkm_filename):
    rpkm_h5 = tables.openFile(rpkm_filename, mode='r')
    d = rpkm_h5.root.rpkm.rpkm.read(field="rpkm")
    rpkm_h5.close()
    return d

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("rpkm_directory", help="Directory of RPKM files to process.")
    parser.add_argument("--negative-binomial", "--nb", help="Use negative binomial model.")
    parser.add_argument("--outfile")
    args = parser.parse_args()

    rpkm_files = filter(lambda x: x.endswith(".h5"), os.listdir(args.rpkm_directory))
    probes = xrange(len(loadRPKM(rpkm_files[0])))
    out_data = pandas.DataFrame(columns=rpkm_files, index=probes)

    # load all the data
    for rpkm_filename in rpkm_files:
        out_data[rpkm_filename] = loadRPKM(rpkm_filename)

    cPickle.dump(out_data, open(args.outfile, 'w'))
    #if args.negative_binomial:
        # apply glm here

    #else:
        # apply zrpkm here

    # apply SVD
    # U, S, Vt = np.linalg.svd(rpkm,full_matrices=False)
    # new_S = np.diag(np.hstack([np.zeros([components_removed]), S[components_removed:]]))
    # # reconstruct data matrix
    # rpkm = np.dot(U, np.dot(new_S, Vt))