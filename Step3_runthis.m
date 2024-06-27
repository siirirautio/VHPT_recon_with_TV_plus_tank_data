%===================================================================================================
%===================================================================================================
% This script computes the Total Variation reconstruction of the deblurred
% CGO sinogram.

% The script calls to tomo_tv.m, a modified version of the primal-dual 
% algorithm with total generalized variation (TGV) tomography.
% copyright (c) 2012 Kristian Bredies

% Copyright (c) 2024 Samuli Siltanen, Siiri Rautio
% Date last modified: June 2024
%===================================================================================================
%===================================================================================================

clear all;
close all;

%===================================================================================================
%========================================== Preliminaries ==========================================
%===================================================================================================

% Add directory with stored subroutines and inputs
addpath Step3_subroutines/

% Load the measurement matrix and parameters from file
load RadonMatrix A measang N P Nang

% Load the deblurred sinogram from file and pad it
filename = 'deblurred_sinogram_0000';
sinogram = double(load(['output/deblurred_sinograms/',filename,'.mat']).deblurred_sinogram);
mn = [sinogram;zeros(1,100)];

% Plot
figure(2)
clf
imagesc(sinogram)
title('Sinogram')
colormap gray
axis off

% Normalize measurement matrix so that its norm is one. Normalize 
% measurement similarly
normA = normest(A);
A = A/normA; 
m = mn/normA; 
m = m/norm(m);

%===================================================================================================
%================= Filtered backprojection reconstruction (for comparison) =========================
%===================================================================================================

% Compute FBP reconstruction
reco_fbp = iradon(m, measang);
reco_fbp = rot90(reco_fbp);

% Plot
figure(3)
clf
imagesc(reco_fbp)
title('FBP reconstruction')
colormap gray
axis off
axis equal

%===================================================================================================
%============================= Total Generalized Variation reconstruction ==========================
%===================================================================================================

% Maximum number of iterations
max_iter = 5000;

% Loop over a range of regularization parameters
for alpha = [0.002, 0.004, 0.006, 0.008]

    % Reconstruct using total generalized variation
    recn = tomo_tv(m, A, 2, alpha, max_iter, 1);
    recn = rot90(recn);

    % From mu to conductivity
    recn_cond = (1-recn)./(1+recn);

    % Plot
    figure()
    clf
    imagesc(reshape(recn_cond,[N,N]))
    title (alpha)
    axis image
    colormap gray
    axis off
    
    % Save as .png file
    reco_plot = uint8(rescale(recn_cond,1,256)); % Rescale for png image
    path = ['output/reconstructions/rec_tv_', filename, 'alpha', num2str(alpha)];
    imwrite(reco_plot,path,'png');
    
    % Save as .mat file
    save(path, 'recn_cond');
end









