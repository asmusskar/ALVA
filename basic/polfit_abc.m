function abc_zr = polfit_abc(E,nu,zi,Lam1,alva)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function evaluates the a, b and c-value for the polynomial in the
% approximation method used (see [1]) to speed up the integration. The 
% output abc-matrix contains the abc constants in two rows. Row 1 the 
% z-related values and row 2 the r-related values

% The polynomial approximations: 
% uz: f1(m) = (a*m^2+b*m+c)*exp(-2*lam1*m)-(2-2*Nu)  % for uz deformations
% ur: f2(m) = (a*m^2+b*m+c2)*exp(-2*lam1*m)-(1-2*Nu) % for ur deformations

% INPUT PARAMETERS:
% Lam1         : relative depth of top layer (=z1/H)
% A, B, C, D   : integral parameters
% xip_z, xip_r : integration points
% nu           : poisson's ratio
% E            : Young's moduli

% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
%[1] Andersen, S., Levenberg, E., & Andersen, M. B (2020). Efficient 
%    reevaluation of surface displacements in a layered elastic half-space. 
%    The International Journal of Pavement Engineering 21(4), 1-8. 
%    https://doi.org/10.1080/10298436.2018.1483502
%--------------------------------------------------------------------------

alpha = [0.01 0.05 0.1]/100;
mz = 1/(2*Lam1)*log(-1./(alpha*(2*nu(1)-2)));
mr = 1/(2*Lam1)*log(-1./(alpha*(2*nu(1)-1)));

% Evaluate A, B, C and D values used for evaulation of a, b and c
ABCD = arb_func_plain(length(E),mz,zi,E,nu,alva);

Azi  = ABCD(1,:)'; % Layer 1 A value
Bzi  = ABCD(2,:)'; %     ... B ...
Czi  = ABCD(3,:)'; %     ... C ...
Dzi  = ABCD(4,:)'; %     ... D ...

ABCD = arb_func_plain(length(E),mr,zi,E,nu,alva);
Ari  = ABCD(1,:)'; % Layer 1 A value
Bri  = ABCD(2,:)'; %     ... B ...
Cri  = ABCD(3,:)'; %     ... C ...
Dri  = ABCD(4,:)'; %     ... D ...

%%%%%%% Evaluate coefficients for polynomial approximations %%%%%%%%
m1  =  mz(1);
m2  =  mz(2);
m3  =  mz(3);

% Denominators
den1 = (m1-m3)*(m1-m2); % denominator 1
den2 = (m2-m3)*(m1-m2); % denominator 2
den3 = (m2-m3)*(m1-m3); % denominator 3

% Organize inverted matrix used in the evaluation abc = Mi*Rhs
Mz = [     exp(2*Lam1*m1)/den1     ,    -exp(2*Lam1*m2)/den2     ,      exp(2*Lam1*m3)/den3
      -exp(2*Lam1*m1)*(m2+m3)/den1 , exp(2*Lam1*m2)*(m1+m3)/den2 , -exp(2*Lam1*m3)*(m1+m2)/den3
        exp(2*Lam1*m1)*m2*m3/den1  , -exp(2*Lam1*m2)*m1*m3/den2  ,   exp(2*Lam1*m3)*m1*m2/den3 ];

% For the radial displacements
m1  =  mr(1);
m2  =  mr(2);
m3  =  mr(3);

% Denominators
den1 = (m1-m3)*(m1-m2); % denominator 1
den2 = (m2-m3)*(m1-m2); % denominator 2
den3 = (m2-m3)*(m1-m3); % denominator 3

% Organize inverted matrix udsed in abc = Mi*Rhs
Mr = [     exp(2*Lam1*m1)/den1     ,    -exp(2*Lam1*m2)/den2     ,      exp(2*Lam1*m3)/den3
      -exp(2*Lam1*m1)*(m2+m3)/den1 , exp(2*Lam1*m2)*(m1+m3)/den2 , -exp(2*Lam1*m3)*(m1+m2)/den3
        exp(2*Lam1*m1)*m2*m3/den1  , -exp(2*Lam1*m2)*m1*m3/den2  ,   exp(2*Lam1*m3)*m1*m2/den3 ];
   
% Evaluate a, b and c coefficients for uz integration
abc_zr = zeros(2,3);

% Evaluate the uz-related paramters
Rhs_z       = (Azi - Czi*(2-4*nu(1))).*exp(-mz'*Lam1) - (Bzi+Dzi*(2-4*nu(1))) + (2-2*nu(1));
abc_zr(1,:) = Mz*Rhs_z;

% Evaluate the ur-related paramters
Rhs_r       = (Ari + Cri).*exp(-mr'*Lam1) + Bri-Dri + (1-2*nu(1));
abc_zr(2,:) = Mr*Rhs_r;

% Set a, b, c equal to zero if E1 approx. E2 approx. = E3*...
if (max(E)-min(E))/max(E)*100 < 5
   abc_zr = zeros(2,3);
end

% *Note: If all layer moduli are equal, i.e., E1 = E2 = E3 = ... we have a
% half-space (one layer) solution and a = b = c = 0. Numerical noise can 
% be experienced if E1 approx E2 approx E3 ... A tolerence is therfore
% introduced, set to 5%, i.e., if 0.95*E1 <= E2 <= 1.05*E1, a, b og c is 0, 
% This is because its mathematically may lead to some isssues determining 
% a, b og c, and because we are so close to the half-space solution. This 
% special case is only used when all layer moduli are approx. equal. If one
% layer moduli differ more than 5%, its treated as a multi-layered model
