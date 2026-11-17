# NAME: convert_to_h5ad.py
# HH
"""
Executes src/01_preprocess/h5ad_utils.py to convert raw data/broad/Sample_<id>.tar(.gz) 
archives to data/converted_h5ad_data/Sample_<id>.h5ad.

"""

from unicodedata import name

from h5ad_utils import convert, load_h5ad_sample, load_h5ad_samples, h5ad_path
# Convert all samples in metadata.json to .h5ad files (skip existing unless overwrite=True)
# Save the .h5ad files in data/converted_h5ad_data/

if name == "__main__":
    convert(overwrite=False)  # Set overwrite=True to rewrite existing .h5ad files
    # step one: set h5ad_path?
    # step two: load all h5ad samples
    # step three: convert to h5ad format