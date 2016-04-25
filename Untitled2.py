
# coding: utf-8

# # Intro
# 
# Name:  
# 
#     star_codes
# 
# Purpose:  
# 
#     To compile all the 4STAR analysis codes within a single module
#     Should be compared to matlab codes created by Yohei
#     contains multiple codes including:
#         - allstarmat - compile raw ascii files into a matlab binary file
#         - starsun = run through typical analysis of direct beam analysis
#         - starzen - for zenith radiance analysis
#         - others...
# 
# Calling Sequence:
# 
#     import star_codes as sc
# 
# Input:
# 
#     See each individual method within the module
# 
# Output:
#    
#     Save matlab file of raw or processed files
#     Some ploting utilities
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - hdf5storage : for saving and reading
#     - Run_Libradtran
#     - os
#     - sys
#     - matplotlib
#     - argparse
#     - Tkinter (for file open dialogs)
#     - load_modis for mat2py_time and toutc
#   
# Needed Files:
# 
#   - raw ascii files or yyyymmddstar.mat file
#   
# Example:
# 
#    to come...
# 
# History:
# 
#     Written: Samuel LeBlanc, Flight from Honolulu to SFO over the Pacific, 2016-04-22

# # Prepare the initial setups

# In[ ]:



