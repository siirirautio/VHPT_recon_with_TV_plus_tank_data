%===================================================================================================
%===================================================================================================
% This function follows the theory in Astala & Päivärinta, 2006 to solve a boundary integral 
% equation computed from DN maps, to compute the boundary traces of CGO solutions 
% omega^+ and omega^-, and construct the integral of their difference over the domain boundary.
% The output is used in the VHPT imaging chain.

% We assume the spatial domain is the unit disc. The mu-Hilbert matrices computed from DN maps 
% are accepted as inputs.

% INPUTS:
% Nang = Number of virtual X-ray angles phi will will compute the solution at
% Nth  = Number of spatial domain boundary points used in solving the BIE
% phi  = Vector of virtual X-ray angles
% tau  = 1D spectral grid corresponding to each angle phi
% Hmu  = Matrix for mu-Hilbert transform
% Hmmu = Matrix for (-mu)-Hilbert transform
% H0   = Matrix for 0-Hilbert transform

% OUTPUTS:
% tildeT_odd = a complex contour integral over the boundary of the difference (omega^+) - (omega^-)

% Authors: Melody Alsaker, Siiri Rautio
% Date last modified: June 2024
%===================================================================================================
%===================================================================================================
function tildeT_odd = compute_CGOs(Nang, Nth, phi, tau, Hmu, Hmmu, H0)

% Complex-valued spatial grid points on the unit disc boundary, used in solving the BIE
th       = 2*pi*(0:(Nth-1))'/Nth;  
zb       = exp(1i*th).'; 

% Extract grid sizes
Nth = numel(th);
Ntau = numel(tau);
Ntrig = (size(Hmu,1)-1)/2;

% Loop over virtual X-ray angles, solving the BIE twice for each angle to get omega^+ and omega^-
parfor aaa = 1:Nang
    kvec = exp(1i*phi(aaa))*tau; % Vertical vector of complex spectral k-values at current angle phi

    % Initialize matrices
    tracemat  = zeros(Nth,Ntau);
    tracematm = zeros(Nth,Ntau);
    wtrace_plus = zeros(Nth,Ntau);
    wtrace_minus = zeros(Nth,Ntau);

    % Computation for +\mu
    for kk = 1:Ntau
        % SOLVE INTEGRAL EQUATION
        % We construct and invert the matrix of operator I - Pmu - P^k_0 in the trigonometric basis
        k = kvec(kk);
        % Initialize operator matrix
        opermat = zeros(2*(2*Ntrig+1));

        % Construct matrix column by column
        for iii = 1:2*(2*Ntrig+1)
            coef           = zeros(2*(2*Ntrig+1),1);
            coef(iii)      = 1;
            Pmuterm        = coef2fun(Pmu_matrix(Ntrig,Hmu,Hmmu)*coef,Ntrig,Nth);
            coef2          = fun2coef(exp(-1i*k*exp(1i*th)).*coef2fun(coef,Ntrig,Nth),Nth,Ntrig);
            Pk0term        = exp(1i*k*exp(1i*th)).*coef2fun(P0_matrix(Ntrig,H0)*coef2,Ntrig,Nth);
            opermat(:,iii) = coef - fun2coef(Pmuterm+Pk0term,Nth,Ntrig);
        end

        expcoef = fun2coef(-exp(1i*k*exp(1i*th)),Nth,Ntrig);
        solcoef = pinv(opermat)*expcoef;
        sol     = coef2fun(solcoef,Ntrig,Nth);

        % Record current result
        tracemat(:,kk) = sol(:);
        w_plus = sol.*exp(-1i*k*exp(1i*th))-1;

        wtrace_plus(:,kk) = w_plus(:);

    end

    % Computation for -\mu
    for kk = 1:Ntau
        % SOLVE INTEGRAL EQUATION
        % We construct and invert the matrix of operator I - Pmu - P^k_0 in the trigonometric basis
        k = kvec(kk);
        % Initialize operator matrix
        opermat = zeros(2*(2*Ntrig+1));
        % Construct matrix column by column
        for iii = 1:2*(2*Ntrig+1)
            coef           = zeros(2*(2*Ntrig+1),1);
            coef(iii)      = 1;
            Pmmuterm       = coef2fun(Pmmu_matrix(Ntrig,Hmu,Hmmu)*coef,Ntrig,Nth);
            coef2          = fun2coef(exp(-1i*k*exp(1i*th)).*coef2fun(coef,Ntrig,Nth),Nth,Ntrig);
            Pk0term        = exp(1i*k*exp(1i*th)).*coef2fun(P0_matrix(Ntrig,H0)*coef2,Ntrig,Nth);
            opermat(:,iii) = coef - fun2coef(Pmmuterm+Pk0term,Nth,Ntrig);
        end

        expcoef = fun2coef(-exp(1i*k*exp(1i*th)),Nth,Ntrig);
        solcoef = pinv(opermat)*expcoef;
        sol     = coef2fun(solcoef,Ntrig,Nth);

        % Record current result
        tracematm(:,kk) = sol(:);
        w_minus = sol.*exp(-1i*k*exp(1i*th))-1;
        wtrace_minus(:,kk) = w_minus(:);

    end

    omega_mat_plus = flip(wtrace_plus.');                   
    omega_mat_minus = flip(wtrace_minus.');            
    omega_diff = omega_mat_plus - omega_mat_minus; 

    % Integrate over the spatial z variable using a complex contour integral over the unit disc
    tildeT_odd(:,aaa) = 1i*(zb*(omega_diff.')/Nth).';

end

tildeT_odd=conj(tildeT_odd)/2; 

end

