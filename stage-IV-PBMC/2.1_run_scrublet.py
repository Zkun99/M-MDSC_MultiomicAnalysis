# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 19:48:28 2024

@author: wllab_bioinformatics
"""

import sys
import os
import scrublet as scr
import scipy.io
import pandas as pd
import matplotlib.pyplot as plt

main_dir = r'E:\NSCLC_PBMC_scRNA_seq\scrublet'

for sub_dir in os.listdir(main_dir):
    sub_dir_path = os.path.join(main_dir, sub_dir)
    
    if os.path.isdir(sub_dir_path):
        print(f"Processing folder: {sub_dir_path}")
        
        counts_matrix = scipy.io.mmread(os.path.join(sub_dir_path, 'matrix.mtx.gz')).T.tocsc()
        out_df = pd.read_csv(os.path.join(sub_dir_path, 'barcodes.tsv.gz'), header=None, index_col=None, names=['barcode'])
        
        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
        
        print(f"Detected doublet rate for {sub_dir}: {scrub.detected_doublet_rate_}")
        
        out_df['doublet_scores'] = doublet_scores
        out_df['predicted_doublets'] = predicted_doublets
        output_file = os.path.join(sub_dir_path, 'doublet.txt')
        out_df.to_csv(output_file, index=False, header=True)
        print(f"Results saved to: {output_file}\n")
        
        plt.figure()
        scrub.plot_histogram()
        hist_file = os.path.join(sub_dir_path, 'doublet_scores_histogram.png')
        plt.savefig(hist_file)
        plt.close()
        print(f"Histogram saved to: {hist_file}\n")

print("All folders processed.")
