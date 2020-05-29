function alva = init_LET(alva)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function initialize the layered elastic half space response
% calculations utilizing Linear Elastic Theory (LET) [1],[2],[3]. The 
% string 'Full' should be used if displacements and stresses inside the 
% medium is considered. For fast surface deformation algorithms utilizing
% the method proposed in [4] the string 'PolFit' should be applied, i.e., 
% if deformations for z = 0 is sought only.

% INPUT PARAMETERS
N  = alva.N;    % Number of Bessel roots in integration
n  = alva.n;    % Number of integration points in between each Bessel root
Xd = alva.Xd;   % Output coordinates
Xl = alva.Xl;   % Load posistion coordinates
a  = alva.a;    % Load radii
E  = alva.E;    % Layer Young's moduli
nu = alva.nu;   % Layer Poisson's ratios
zi = alva.zi;   % Layer interface depths (i.e. depth of layer n-1 out
                % of n layers (layer n depth goes to infinity)
% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
%[1] Burmister, D. M. (1945). The general theory of stresses and 
%    displacements in layered systems. i. Journal of applied physics, 16(2)
%    , 89–94.
%[2] Huang, Y.H., 2003. Pavement Analysis Design, 2nd Edition, 
%    Prentice-Hall, New Jersey.
%[3] Ioannides, A. M., & Khazanovich, L. (1998). General formulation for 
%    multilayered pavement systems. Journal of transportation engineering, 
%    124(1), 82–90.
%[4] Andersen, S., Levenberg, E., & Andersen, M. B (2020). Efficient 
%    reevaluation of surface displacements in a layered elastic half-space. 
%    The International Journal of Pavement Engineering 21(4), 1-8. 
%    https://doi.org/10.1080/10298436.2018.1483502
%--------------------------------------------------------------------------

% General parameters
H    = zi(end);       alva.H    = H;     % Bottom depth of last finite layer
Lam1 = zi(1)/zi(end); alva.Lam1 = Lam1;  % Relative height of top layer

% Integration points and weights
XipWip = numint_coeff(N,n,Xd,Xl,a,H); alva.XipWip = XipWip;

if strcmp(alva.analysis,'Full')       % Full response analysis
    
    % Evaluate coefficients of integration
    alva.ABCD = arb_func(length(E),zi,E,nu,alva);
 
    % Evaluate response 
    alva = LET_response(alva);  

elseif strcmp(alva.analysis,'PolFit') % Approx. polynomial fit method
    if sum(Xd(:,3)) > 0
        error('Change output locations to surface points')
    end
    
    % Integration points for the polynomial integrals are. These are high.
    N0      = 1000; 
    n0      = 30;
    XipWip0 = numint_coeff_plain(N0,n0,Xd,Xl,a,H);
    
    % Evaluate coefficient proportional integrals
    alva.I = polfit_int(Xd,Xl,a,H,Lam1,XipWip0,nu(1));
        
    % Evaluate surface displacements
    alva   = LET_response_polfit(alva);

else
    error('Analysis type unknown')
end

