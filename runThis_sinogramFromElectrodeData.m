%===================================================================================================
%===================================================================================================
% This script loads EIT data, computes a DN map, and calls subroutines necessary to generate a
% blurred CGO sinogram from real-world EIT electrode data, as part of the VHPT imaging chain.
% Output from this routine serves as an input into a neural network which strips the higher-order
% terms of the scattering series and peforms deconvolution to achieve a sharpened sinogram.
%
% This code has not been optimized for speed. 
%
% If plotting options are specified, we will compute a simple FBP reconstruction from the blurred
% sinogram, to provide a baseline check that the algorithm is working properly.
%
% In this routine we assume:
% -The spatial domain is the unit disc
% -EIT measurements are taken using an EVEN number of electrodes L
% -Current is applied on L electrodes according to a trigonometric current pattern basis as used
%  by the ACT5 EIT machine (Shishvan, Abdelwahab, et al, 2023).

% To calibrate the DN map, we use the approach based on the work of Garde-Hyvönen, 2021, using
% voltages collected on a homogeneous saline tank as well as data simulated on a homogeneous domain,
% computed using the Complete Electrode Model (CEM) solved via FEM. These three datasets must be
% specified and loaded.

% Authors: Melody Alsaker, Siiri Rautio
% Date last modified: June 2024
%===================================================================================================
%===================================================================================================

clear 
totalRuntimeStart = tic;    % Start code timer

%===================================================================================================
%====================================== User-Specified Options =====================================
%===================================================================================================

% Do you wish to use CPU parallelization over the VHPT angles?
parallelize = true;

% Specify options for saving and plotting. We'll plot the sinogram and a simple FBP reconstruction.
save_sinogram           = true;  % Save the CGO singogram to a .mat file?
plot_to_screen          = true;  % Display a plot of the CGO sinogram and FBP recon to screen?
save_plots              = true;  % Save a plot of the CGO sinogram and FBP recon as a .png file?

% Specifiy directory where output will be saved. Include forward or backslash as appropriate.
% If it doesn't exist already, we'll create it
outdir  = 'output/CGO_sinograms/';

% Specify directory where all EIT data is stored
datadir = 'voltage_data/';

% Specify names of .mat files containing EIT data
V_clb_fname = 'saline_270_22_11_07_10_16_18';        % Calibration data (e.g. homogeneous tank)
V_trg_fname = 'two_large_yellow_22_11_07_10_21_29';  % Target data (what we wish to image)
V_1cem_fname = 'CEM_voltages_homog_CPamp0p35';       % Simulated homogeneous data from CEM


% Specify computational parameters
tauMAX      = 6;        % Truncation radius for spectral frequency parameter (k-grid)
Ntau_half   = 20;       % Spectral grid size: we'll have 2x this many tau values for each angle phi
Nang        = 100;      % Number of VHPT angles phi in [0, 2*pi) (number of columns in sinogram)
ee          = 0.05;     % Used to compute width of Gaussian window in FT. Smaller = more truncation
t_pseudoMAX = 3;        % pseudo-time range will be [-t_pseudoMAX, t_pseudoMAX]
Nt_pseudo   = 256;      % Number of points in 1D pseudotime grid (number of rows of sinogram)
Nth         = 128;      % Number of spatial domain boundary points used in solving the BIE


%===================================================================================================
%========================================== Construct DN Maps ======================================
%===================================================================================================
[DN, DN1, Ntrig] = compute_DN_maps(datadir, V_clb_fname, V_trg_fname, V_1cem_fname);


%===================================================================================================
%================================ Compute CGO Traces over Disc Boundary ============================
%===================================================================================================

% Construct mu-Hilbert matrices used in construction of projection operators
[H0,Hmu,Hmmu] = muHilbert(DN,DN1,Ntrig);

% Get weights and abscissae for Gaussian quadrature over 1D spectral tau grid
[tau_neg,tauWeights_neg] = gaussint(Ntau_half,-tauMAX,0);
[tau_pos,tauWeights_pos] = gaussint(Ntau_half,0,tauMAX);
tau  = [tau_neg(:);tau_pos(:)];
tauWeights = [tauWeights_neg(:);tauWeights_pos(:)];

% Define the set of Nang virtual X-ray measurement angles phi in [0, 2pi)
phi = (0:(Nang-1))/(Nang-1)*2*pi;

% If parallelization is desired, open parallel pool
if  parallelize == true
    if isempty(gcp('nocreate'))
        parpool
    end
end

% Solve the BIE to compute the boundary traces of CGO solutions omega^+ and omega^-
% and construct the integral of their difference over the domain boundary.
fprintf('Computing the CGO traces for %u virtual X-ray angles phi... \n \n', Nang);
tildeT_odd = compute_CGOs(Nang, Nth, phi, tau.', Hmu, Hmmu, H0);


%===================================================================================================
%================================ Apply Windowed Fourier Transform  ================================
%===================================================================================================

% Columnwise application of a 1D Fourier Transform transforms from the spectral frequency parameter 
% tau into the pseudo-time parameter t. We use a windowing to truncate noisy spectral data in the 
% high spectral frequencies which stems from ill-posedness

% Create 1D pseudo-time grid
t_pseudo  = linspace(-t_pseudoMAX,t_pseudoMAX,Nt_pseudo);
Dt_pseudo = t_pseudo(2)-t_pseudo(1);

% Construct Gaussian window function automatically adjusted to tauMAX
a       = -log(ee)/tauMAX^2;
g_window = exp(-a*tau.^2);

% Initialize matrix
T_odd    = zeros(Nt_pseudo,Nang);


% Loop over X-ray directions phi, taking 1D windowed Fourier transform of tildeT_odd. 
% This transforms from the spectral domain to the pseudo-time domain
for phii = 1:Nang
    for iii = 1:Nt_pseudo
        T_odd(iii,phii) =...
            (tauWeights.') * (exp(-1i*t_pseudo(iii)*tau) .* g_window.*tildeT_odd(:,phii));
    end
end

%===================================================================================================
%========================= Compute Final Sinogram and FBP Reconstruction  ==========================
%===================================================================================================

% Apply complex phase function
explfun = repmat(-2*pi*1i*exp(1i*phi), Nt_pseudo,1);
T_odd_exp    = T_odd./explfun;

% Apply fundamental theorem of calculus to undo derivative, via columnwise integration
sinogram = Dt_pseudo * cumsum(T_odd_exp,1);
sinogram=flip(real(sinogram)); % Final sinogram!

% Compute a simple FBP reconstruction of the blurry sinogram as a quick check
mu = iradon(sinogram, phi/(2*pi)*360);  % Beltrami paramater 
sigma = (1-mu)./(1+mu);

totalRuntime = toc(totalRuntimeStart);
fprintf('Total program runtime: %f seconds \n \n',totalRuntime);


%===================================================================================================
%=================================== Saving and Plotting Routines  =================================
%===================================================================================================

if save_sinogram == true
    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    outFname = [outdir, V_trg_fname, '_tauMAX', strrep(num2str(tauMAX),'.','p'), ...
                                     '_Nang', num2str(Nang), ...
                                     '_tpseudoMAX', strrep(num2str(t_pseudoMAX),'.','p'),  ];
    save(outFname, 'sinogram');
end


if (plot_to_screen == true || save_plots == true)

    if( plot_to_screen == true )
        figure
    else
        figure('visible','off')
    end

    t = tiledlayout(1,2, 'TileSpacing','Compact','Padding','Compact');

    n1 = nexttile;
    imagesc(sinogram)
    axis square
    set(gca,'XTick',[]); set(gca,'YTick',[])
    title('CGO Sinogram $R_{\mathrm{odd}}\mu(s,\varphi)$', 'interpreter', 'latex')
    colormap(n1, gray)

    n2 = nexttile;
    imagesc(flipud(sigma))
    axis square
    set(gca,'XTick',[]); set(gca,'YTick',[])
    title('Simple \sigma recon from FBP')
    colormap(n2, gray)

    saveas(gcf,[outFname,'.png'])

end


