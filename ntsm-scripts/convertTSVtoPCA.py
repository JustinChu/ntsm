import argparse
from pyfaidx import Fasta
import os
import numpy as np
import sys
from sklearn import decomposition
from sklearn import datasets
from sklearn.decomposition import PCA
import pandas as pd

class GenPCA:
    def __init__(self, matrix, numComp):
        self._matrix = matrix
        self._n = numComp
    
    def runPCA(self, prefix):
        dataset = pd.read_csv(self._matrix, sep='\t', header = 0, index_col = 0)
        pca = PCA(n_components=self._n)
        X_reduced = pca.fit_transform(dataset.transpose())
        df = pd.DataFrame(pca.components_.transpose())
        df = df.set_axis(pca.feature_names_in_, axis=0)
        df.to_csv(prefix + "_rotationalMatrix.tsv", sep='\t', index_label = "AlleleID")
        pc = pd.DataFrame(X_reduced)
        pc = pc.set_axis(dataset.axes[1], axis=0)
        pc.to_csv(prefix + "_components.tsv", sep='\t', index_label = "SampleID")


def main():
    parser = argparse.ArgumentParser(description='Computes a PCA on various matrix produced by tools like ntsmVCF. Produces a rotational matrix for n components.')
    parser.add_argument("-m", '--matrix', type=str, dest='matrix', 
                        help='input matrix', required=True)
    parser.add_argument("-n", '--numComp', type=str, dest='numComp', 
                        help='Number of comp', required=False, default = 20)
    parser.add_argument("-p", '--prefix', type=str, dest='prefix', 
                        help='Output prefix', default = "")

    args = parser.parse_args()
    gen = GenPCA(args.matrix, args.numComp)
    gen.runPCA(args.prefix)

main()