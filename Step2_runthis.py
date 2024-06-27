#========================================================================================================
#========================================================================================================
# This script loads the blurry CGO sinograms and neural network model parameters and then deblurs the 
# sinograms while simultaneously stripping away the higher order terms from the scattering series. 
# The blurry sinograms are computed using the 'Step1_runthis' routine. 
# The output is saved as a .mat file.

# Copyright (c) 2024 Siiri Rautio
# Last modified: April 2024
#========================================================================================================
#========================================================================================================

#========================================================================================================
#===================================== Import packages ==================================================
#========================================================================================================

import cv2
import glob
import scipy.io
from tensorflow.keras import models
import numpy as np
from scipy.io import savemat

#========================================================================================================
#========================================== Preliminaries ===============================================
#========================================================================================================

# Data dimensions
M = 100
N = 200
X = 1

test_input = []

#========================================================================================================
#======================================= Load and preprocess data =======================================
#========================================================================================================

# Path to files
files      = sorted(glob.glob ("output/CGO_sinograms/*.mat"))

for myFile in files:
    print(myFile)
    mat    = scipy.io.loadmat(myFile)
    matrix = mat['sinogram']
    
    # Normalize data
    min_input = np.amin(matrix)
    max_input = np.amax(matrix)
    matrix    = (matrix - min_input) / (max_input - min_input)
    
    # Resize data
    test_input.append (cv2.resize(matrix.real, dsize=(M,N), interpolation=cv2.INTER_CUBIC))
test_input = np.array(test_input)    

# Reshape data
test_input = test_input.reshape(-1, N, M, X)
test_input = test_input.astype('float64')

#========================================================================================================
#=========================================== Deblur sinograms ===========================================
#========================================================================================================

# Load model 
autoencoder = models.load_model('Step3_models/model_140224')

# Neural network prediction
test_prediction = autoencoder.predict(test_input)

# Save results as mat.-files
for iii in range(len(test_prediction)):
    myData   = test_prediction[iii,:,:,:]
    index    = iii + 1
    filename = "output/deblurred_sinograms/S%i_deblurred_sinogram.mat" % index
    savemat(filename, {'deblurred_sinogram': myData})
    print(index)

