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

### Data

All the data for the 3 phantoms used in the paper are included:
- pac_man_2_22_11_07_11_11_05.mat
- split_yellow_22_11_07_10_33_25.mat
- two_large_yellow_22_11_07_10_21_29.mat

### How to run the code

Run the file `runThis_sinogramFromElectrodeData.m`. The first time you run the code it will generate an output directory if the option to save is selected. This file outputs the blurry sinogram, that is the input for the neural network.

Run the file `deblur_sinogram.py` to deblur the sinogram.

Finally, run the file `VHPT_TV_comp.m` to reconstruct the deblurred sinogram using Total Variation regularization and compute the reconstucted conductivity.
