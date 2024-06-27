# Computed tomography without X-rays: parallel-beam imaging from nonlinear current flows

### General information

This code was used to produce the results in "Computed tomography without X-rays: parallel-beam imaging from nonlinear current flows" by Melody Alsaker*<sup>1</sup>, Siiri Rautio<sup>2</sup>, Fernando Moura<sup>3</sup>,  Juan Pablo Agnelli<sup>4</sup>, Rashmi Murthy<sup>5</sup>, Matti Lassas<sup>2</sup>, Jennifer L. Mueller<sup>6</sup>, Samuli Siltanen<sup>2</sup>.

1. Department of Mathematics, Gonzaga University, Spokane, WA 99258, USA
2. Department of Mathematics and Statistics, University of Helsinki, Helsinki, Finland
3. Engineering, Modeling and Applied Social Sciences Center, Federal University of ABC, São Paulo, Brazil
4. FaMAF, National University of Córdoba, Córdoba, Argentina and CIEM, National Scientific and Technical Research Council (CONICET), Argentina
5. Department of Mathematics, Bangalore University, India
6. Department of Mathematics & School of Biomedical Engineering, Colorado State University, Fort Collins, CO 80521, USA

*Corresponding author. Email:  alsaker@gonzaga.edu.

### Code authors

Main algorithms: copyright (c) 2004 Melody Alsaker, Siiri Rautio, Samuli Siltanen

Some subroutines use code authored by other individuals. Whenever possible, appropriate credit has been given within the relevant files. 

### Data

The EIT data for the 3 targets used in the paper, along with the calibration set (collected on a tank filled with homogeneous saline) and the set corresponding to numerically simulated data on a homogeneous domain are included in the directory labeled DATA:
- Data_S1: Target I 
- Data_S2: Target II 
- Data_S3: Target III 
- Data_S4_clb: Calibration data (homogeneous tank)
- Data_S5_cem: CEM simulated data (homogeneous domain)

### How to run the code

1. Run the MATLAB file `Step1_runthis.m`. The first time you run the code it will generate an output directory if the option to save is selected. This file reads in EIT data and outputs a .mat file containing the blurry sinogram, which is the input for the neural network in Step 2.

2. Run the Python file `Step2_runthis.py` to deblur the sinogram which was output in Step 1. The deblurred sinogram is saved as a .mat file in the output directory. 

3. Finally, run the file `Step3_runthis.m` to reconstruct a conductivity image from the deblurred sinogram using Total Variation regularization. The final reconstruction is saved 
