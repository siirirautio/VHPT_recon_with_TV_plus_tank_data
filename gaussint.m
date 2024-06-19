% FUNCTION program [X,W] = GAUSSINT(n,itype,a,b,alfa,beta)
% ----------------------------------------------------------------------
% Function to evaluate Gaussian integration weights and abscissae
% INPUT  n        : number of division points (= order of orthog poly)
%        itype    : type of polynomial and weighting function to be used
%        a,b      : optional interval limits (for itype = 1,...,4)
%        alfa,beta: powers for Jacobi weighting function (itype = 3)
% OUTPUT X(n)     : column vector of Gaussian abscissae (incr order)
%        W(n)     : corresponding column vector of weights
% Itype  Polynomial                   w(x)                default [a,b]
%   1    Legendre  P(n,x)             1                      [-1,1]
%   2    Jacobi    P(n,x,1,1)         1-x*x                  [-1,1]
%   3    Jacobi    P(n,x,alfa,beta)   (1-x)^alfa*(1+x)^beta  [-1,1]
%   4    Chebyshev T(n,x)             1/sqrt(1-x*x)          [-1,1]
%   5    Laguerre  L(n,x)             exp(-x)                [0,inf]
%   6    Hermite   H(n,x)             exp(-x*x)              [-inf,inf]
% -----------------------------------------------------------------------
% DEFAULT INPUT : Gaussint(n)     - w(x)=1 over [-1,1]
%                 Gaussint(n,a,b) - w(x)=1 over [a,b]
% CALLS TO      : none
% 14/11/90      : Addie Lambert - Rolf Nevanlinna Institute
% -----------------------------------------------------------------------

function [X,W] = gaussint(n,itype,a,b,alfa,beta)

maxtyp =6;
scaling=0;                                     % Default : no scaling
if (nargin==1)  itype = 1;  end;               % GAUSSINT(n)
   if (nargin==3)                              % GAUSSINT(n,a,b)-scaling
      b=a;  a=itype;  itype=1;  scaling=1;
   elseif (nargin>2) & (itype<5)               % [a,b] specified-scaling
      scaling=1;
   elseif (nargin>2) & (itype>4)
      disp('No interval scaling allowed for ITYPE=5 or 6')
      disp('--- Values of a and b have been ignored')
   end;
   if (itype>maxtyp) | (itype<1) ...
      error('Error in weighting type, gaussint(n,ITYPE,..)');
   end;
%                       --------------------
% Set up gamsq, delta and factor for the selected polynomial type
if (itype==1)                                  % Legendre:P(n,x)
   m     = (1:n);
   gamsq = (m-1).^2 ./ ((2*m-1).*(2*m-3));
   delta = zeros(n,1)';
   factor= 2;  scalef=1;
elseif (itype==2)                              % Jacobi:P(n,x,1,1)
   m     = (1:n);
   gamsq = (m.^2-1) ./ (4*m.^2-1);
   delta = zeros(n,1)';
   factor= 4/3;
   if (scaling==1)  scalef=(b-a)^2/4;  end;
elseif (itype==3)                              % Jacobi:P(n,x,alfa,beta)
   if (nargin<6)
      error('Jacobi:input is GAUSSINT(n,3,a,b,alfa,beta)');
   end;
   m     = (1:n);
   gamsq = 4*(m-1).*(m+alfa-1).*(m+beta-1).*(m+alfa+beta-1);
   denom = (2*m+alfa+beta-1).*(2*m+alfa+beta-2).^2 .*(2*m+alfa+beta-3);
   gamsq = gamsq ./ denom;
   delta = -(alfa^2 - beta^2) ./ ((2*m+alfa+beta-2).*(2*m+alfa+beta));
   factor= 2^(1+alfa+beta) * gamma(alfa+1) * gamma(beta+1) ...
         / gamma(alfa+beta+2);
   if (scaling==1)  scalef=((b-a)/2)^(alfa+beta);  end;
elseif (itype==4)                              % Chebyshev:T(n,x)
   m     = (1:n);                              % -different from other
   X     = cos((2.*m-1) * pi / (2*n))';        % types since X and the
   X     = sort(X);                            % equal weights can be
   W     = ones(n,1) * pi / n;                 % calculated directly
   if (scaling==1) X=((b-a).*X+(b+a))/2; end;
elseif (itype==5)                              % Laguerre:L(n,x)
   m     = (1:n);       factor= 1;
   gamsq = (m-1).^2;    delta = 2*m-1;
elseif (itype==6)                              % Hermite:H(n,x)
   m     = (1:n);       factor= sqrt(pi);
   gamsq = (m-1) ./ 2;  delta = zeros(n,1)';
end
%                       --------------------
% For Chebyshev (itype=4), X and W have already been calculated
% Otherwise, solve eigensystem Jv = Xv to get X and then W
if (itype ~= 4)
   gam  = sqrt(gamsq);
   J    = diag(delta)+diag(gam(2:n),1)+diag(gam(2:n),-1);
   [V,X]= eig(J);
   X    = diag(X);    [X,I] = sort(X);
   W    = zeros(n,1);
   for i = 1:n
      v   = V(:,i) ; v = v/norm(v);
      W(i)= v(1)^2 * factor;
   end;
   W = W(I) ;
   if (scaling==1)                             % If [a,b] specified
      X = ((b-a).*X+b+a) * 0.5;                % - scale X and W
      W = W * scalef * (b-a) * 0.5;            % accordingly
   end;
end;