function alva = LET_response_polfit(alva)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function evaluates the deformations of a layered elastic half space
% model at the surface only, i.e. for z = 0. The approximation method
% proposed by  Andersen et al. (2018) is employed.

% INPUT PARAMETERS:
Xd     = alva.Xd;     % Output coordinates (x,y,z)
Xl     = alva.Xl;     % Load posistion coordinates (x,y)
a      = alva.a;      % Load radii
q      = alva.q;      % Load pressure
E      = alva.E;      % Layer Young's moduli
nu     = alva.nu;     % Layer Poisson's ratios
zi     = alva.zi;     % Layer interface depths (i.e. depth of layer n-1 out
                      % of n layers (layer n depth goes to infinity)
XipWip = alva.XipWip; % Integration points and weights used in integration
I      = alva.I;      % Coefficient proportional integrals
%-------------------------------------------------------------------------

% If zi has more entries than length(E)-1, these are removed. In principle
% the las entry in zi is infinite. However, the infinite value is not
% included in this code, as we do not operate with infinite number. We
% handle this in an optional way.
if length(zi) <= length(E)
    zi(length(E):end) = [];
end

% The code is organized such that we need a minimum of 2 layers! Check if
% zi has at least one value
if isempty(zi)
    display('the code is arranged such that we need minimum two layers!!')
end

% Define depth of layer n-1 out of n layers (layer n depth goes to infinity)
H = zi(end);

% Check if the length of E and nu are equal
if length(E) ~= length(nu)
    display('length(E) ~= length(nu) !!!')
end

% Number of loads (xl) and deformation points (xd)
xl = size(Xl,1); % Number of load points
xd = size(Xd,1); % Number of deformation points

% Make sure that the length of q and a are correct. These should be equal
% to the number of loads.
if length(q) < xl || length(a) < xl
    display('Number of q and/or a values is lower than the number of load coordinates (= number of loads), so all loads are given the same q- and a-values')
    q = q(1)*ones(xl,1);
    a = a(1)*ones(xl,1);
end

% Make sure that q and a are vectors "standing up"
if size(q,1) < size(q,2)
    q = q';
end

if size(a,1) < size(a,2)
    a = a';
end

% Make sure that zi is a 'horizontal' vector
if size(zi,1) > size(zi,2)
    zi = zi';
end

% Increase Xl, q, alpha, Xd and z in order to estimate the radius from each
% of the loads to all of the points. The number of radii are: xl*xd
Xl    = repmat(Xl,[xd,1]);   % Extends with the number of evaluation
                             % points
q     = repmat(q,[xd,1]);    % -||-
alpha = repmat(a/H,[xd,1]);  % -||-
Xd    = repmat(Xd,[xl,1]);   % -||- <-- is reorganized below to 
                             % correspond to the right Xl content/sequence

% Reorganize Xd so the first xl rows correspond to the same deformation
% point, and following xl rows correspond to the next deformation point
% etc... see description of the vector r nine lines below to understand
% the organization
dx = repmat(1:xd,xl,1);
dx = dx(:)' + repmat(0:xd:xd*(xl-1),1,xd);
Xd = Xd(dx,:);

% Evaluate the radii between deformation points and loads
% r = [  deformation point 1 and load 1
%        deformation point 1 and load 2
%        deformation point 1 and load 3
%                        .
%        deformation point 1 and load N
%        deformation point 2 and load 1
%        deformation point 2 and load 2
%                        .
%                        .
%      last deformation point and load N];

r   = sqrt((Xd(:,1)-Xl(:,1)).^2 + (Xd(:,2)-Xl(:,2)).^2);

% Introduce the parameter rho (notice: since r is a vector, rho is too)
rho = r/H;

% Integration points and weights
lrho  = length(rho); % Number of load-to-displacement-points
xip_z = XipWip(0*lrho+1:1*lrho,:);
wip_z = XipWip(1*lrho+1:2*lrho,:);
xip_r = XipWip(2*lrho+1:3*lrho,:);
wip_r = XipWip(3*lrho+1:4*lrho,:);

% Number of columns in reorganized integration points and weights (=number
% of integration points per integral)
nz = size(xip_z,2); % = size(xip_r,2);

% Interpolate A, B, C and D values used for integration
[Az,Bz,Cz,Dz,Ar,Br,Cr,Dr]=arb_func_polfit(E,nu,zi,xip_z,xip_r,lrho,alva);

% Poisson's ratio and Youngs Modulus for the top layer
Nu   = nu(1);
Ee   = E(1);
Lam1 = zi(1)/H;

% Evaluate abc values
abc_zr = polfit_abc(E,nu,zi,Lam1,alva); % ,ABCD added

abc_z = abc_zr(1,:);
abc_r = abc_zr(2,:);

% Organize polynomials for the approximation of the A, B, C and D coefficients
polz  = (abc_z(1).*xip_z.^2+abc_z(2).*xip_z+abc_z(3)).*exp(-2*Lam1*xip_z)-(2-2*Nu);
polr  = (abc_r(1).*xip_r.^2+abc_r(2).*xip_r+abc_r(3)).*exp(-2*Lam1*xip_r)-(1-2*Nu);

% Speeded up solution for the radial displacements, ur:
ur = besselj(1,xip_r.*repmat(rho,1,nz)).*...
    (...
    (Ar+Cr).*exp(-xip_r*Lam1)+...
    (Br-Dr)-polr...
    );

ur  = q.*alpha.*H*(1+Nu)/Ee.*sum(ur.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r,2);

% Evaluate polynimial integral contribution
urp = q/Ee.*(abc_r(1)*I(:,5)+abc_r(2)*I(:,6)+abc_r(3)*I(:,7)+I(:,8));

% Add polynomial integrals to ur
ur = ur + urp;

% Polynomial fit solution
uz = besselj(0,xip_z.*repmat(rho,1,nz)).*...
    (...
    (Az-Cz.*(2-4*Nu)).*exp(-xip_z.*Lam1)-...
    (Bz+Dz.*(2-4*Nu))-polz...
    );

uz =  q.*alpha.*(-1)*H.*(1+Nu)./Ee.*...
    sum(uz.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z,2);

% Evaluate polynomial integral contributions
uzp = q/Ee.*(abc_z(1)*I(:,1)+abc_z(2)*I(:,2)+abc_z(3)*I(:,3)+I(:,4));

% Add polynomial solution to uz
uz = uz + uzp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Transformation between (r,theta) and (x,y) co-ordinates %%%%%%%%%%%%%
% Below a transformation of each of the deformations evaluated above from
% the (r,theta) coordinates to the (x,y) coordinates are made. After the
% transformation, the deformations are added together.

%                         p  (load point)
%   (z is downwards)
%    z-> x              ^ (deformation vector ur; pointing towards the load)
%    |                 /
%  y v                /
%                    o  (deformation point)
%
% NB: the r-axis (goes in the direction


% Transformation matrix (each point-to-load radius results in 3 deformations)
T = sparse(3*length(r),3*length(r),(3*3-4)*length(r));

% Organize content of the T matrix
COS = (Xd(:,1)-Xl(:,1))./r;
SIN = (Xd(:,2)-Xl(:,2))./r;  % Minus on y-coordinate because z-axis points
                             % downward

% Replace content in COS and SIN with, respectively, 1 and 0 if r = 0
COS(find(r==0)) = 1;
SIN(find(r==0)) = 0;

% Transformation matrix is denoted S (and we use its transpose = inverse 
% in the transformation). This is given as:
% S = T =[ cos(theta) sin(theta) 0   = [  Xd(:,1)-Xl(:,1)  ,  Xd(:,2)-Xl(:,2)  ,  0
%         -sin(theta) cos(theta) 0      -(Xd(:,2)-Xl(:,2)) ,  Xd(:,1)-Xl(:,1)  ,  0
%              0           0     1];             0         ,         0         ,  r]*1/r;

% Address the content into the matrix
T(1:3*3*length(r)+3:end)               =  COS;
T(2:3*3*length(r)+3:end)               = -SIN;
T(3*length(r)+1:3*3*length(r)+3:end)   =  SIN;
T(3*length(r)+2:3*3*length(r)+3:end)   =  COS;
T(2*3*length(r)+3:3*3*length(r)+3:end) =  1;

% Organize displacement vector
u1 = zeros(length(r)*3,1);
u1(1:3:end) = ur;
%u1(2:3:end) = 0; % Displacements transverse to the radius (u_theta = 0)
u1(3:3:end) = uz;

% Transform deformations in to (x,y,z)-coordinates
u1=T'*u1;
u=u1;

% Organize matrix Txyz for adding displacements together
unos = zeros(3*xd*xl,1);                % Number of ones in Txyz
unos(1:xl)=1:9*xd:9*xl*xd;
unos(xl+1:2*xl)=unos(1:xl)+3*xd+1;
unos(2*xl+1:3*xl)=unos(1:xl)+2*3*xd+2;

for i=1:xd-1
    unos(3*xl+1+3*xl*(i-1):3*xl+3*xl*i)=unos(1:3*xl)+i*(3*xd*3*xl+3);
end

% Matrix adding contributions together from different loads
Txyz = spalloc(3*xd,3*length(r),length(unos));
Txyz(unos) = 1;

% Final displacement vector
u = Txyz*u;

alva.ux    = u(1:3:end-2);
alva.uy    = u(2:3:end-1);
alva.uz    = u(3:3:end);  