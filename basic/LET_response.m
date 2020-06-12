function alva = LET_response(alva)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function evaluates the response of a layered elastic half-space
% model utilizing Linear Elastic Theory (LET) [1],[2],[3].

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
ABCD   = alva.ABCD;   % Coefficients of integration 

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
% -------------------------------------------------------------------------
% If we only have one layer, we generate two layers

% If zi has more entries than length(E)-1, these are removed. In principle
% the last entry in zi is infinite. However, the infinite value is not
% included in this code, as we do not operate with infinite number. We
% handle this in an optional way.
if length(zi) <= length(E)
    zi(length(E):end) = [];
end

% The code is organized such that we need a minimum of 2 layers! Check if
% zi has at least one value
if isempty(zi)
    disp('the code is arranged such that we need minimum two layers!!')
end

% Define depth of layer n-1 out of n layers (layer n depth goes to infinity)
H = zi(end);

% Check if the length of E and nu are equal
if length(E) ~= length(nu)
    disp('length(E) ~= length(nu) !!!')
end

% Number of loads (xl) and deformation points (xd)
xl = size(Xl,1); % Number of load points
xd = size(Xd,1); % Number of deformation points

% Make sure that the length of q and a are correct. These should be equal
% to the number of loads.
if length(q) < xl || length(a) < xl
    disp('Number of q and/or a values is lower than the number of load coordinates (= number of loads), so all loads are given the same q- and a-values')
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
% etc... see description of the vector r nine lines below to understand the 
% organization

% Built dx vector used to reorganize/sort Xd in the right order
dx = repmat(1:xd,xl,1);
dx = dx(:)' + repmat(0:xd:xd*(xl-1),1,xd);
Xd = Xd(dx,:);
z  = Xd(:,3);

% Evaluate the radii between evaluation point and load
% r = [  evaluation point 1 and load 1
%        evaluation point 1 and load 2
%        evaluation point 1 and load 3
%                        .
%        evaluation point 1 and load N
%        evaluation point 2 and load 1
%        evaluation point 2 and load 2
%                        .
%                        .
%      last evaluation point and load N];

r    = sqrt((Xd(:,1)-Xl(:,1)).^2 + (Xd(:,2)-Xl(:,2)).^2);

% Introduce the parameter rho (notice: since r is a vector, rho is too)
rho = r/H;
rho(rho==0) = 1e-20; % insert a small value to avoid singularity in 
                     % calculation of stresses

% Integration points and weights (takes out the length(rho) rows in XipWip
% related to one load for one position)
lrho = length(rho);                % Number of load-to-displacement-points
xip_z = XipWip(0*lrho+1:1*lrho,:); % takes out the first length(rho) rows
wip_z = XipWip(1*lrho+1:2*lrho,:); % takes out next
xip_r = XipWip(2*lrho+1:3*lrho,:); % etc...
wip_r = XipWip(3*lrho+1:4*lrho,:); % etc...

% Number of columns in reorganized integration points and weights (= number
% of integration points per integral)
nz = size(xip_z,2); 

% The depth z defines which E, nu, Ai, Bi, Ci and Di values should be used
% in the different required evaluations of the response. A
% horizontal vector called layer_no gives the information of which
% values should be used. The length(layer_no) corresponds to the number of
% radii, r, we have (i.e. the length of rho = r/H).

% Layer number where the point we use is given*
layer_no = sum(repmat(z,1,length(zi)) > repmat(zi,length(rho),1),2)'+1; 

% *Note: if (Xd,z) is larger than (zi,1) the function gives 1+1
% = layer no. 2, if (Xd,z) is larger than (zi,2) the function gives 1+1+1 
% = layer no. 3 etc.

Lami  = zeros(length(rho),1);
Lami1 = Lami;
for i = 1:length(rho) % length(layer_no) = n3
    
    % Evaluate the lam_i value (if we are in the bottom layer, Lam_i=infinite)
    if layer_no(i) > length(zi)
        Lami(i) = zi(end)/H*1e20; % In principle Lami is multiplied by 
                                  % A = 0 and C = 0 at the bottom layer, 
                                  % and this value and will not have any 
                                  % contribution
    else
        Lami(i) = zi(layer_no(i))/H;
    end
    
    % Evaluate the lam_{i-1} value (NB! lam_{0} = 0)
    if layer_no(i) == 1 % (= lam_{i-1} = lam{0})
        Lami1(i) = 0;
    else
        Lami1(i) = zi(layer_no(i)-1)/H;
    end
end

% Evaluating response

% Summation form
% length(rho) = number of point to load radii, nz = number of integration 
% points in every point to load integral
Aa = zeros(length(rho),nz); Bb = Aa; Cc = Aa; Dd = Aa;  

for j = 1:length(rho) % number of point to load radii
    % n in function below = number of layers = length(E)
    % m in function below = number of integration points = xip_z(j,:)

    ABCDz = arb_func_interp(length(E),xip_z(j,:),ABCD);
    
    % Select the relevant A, B, C and D values for the response points
    % considered
    Aa(j,:) = ABCDz((layer_no(j)-1)*4+1,:); 
    Bb(j,:) = ABCDz((layer_no(j)-1)*4+2,:); 
    Cc(j,:) = ABCDz((layer_no(j)-1)*4+3,:);
    Dd(j,:) = ABCDz((layer_no(j)-1)*4+4,:);
end

% Define E and Nu values for the layer considered for the evaluation of response
Nu = nu(layer_no)';
Ee = E(layer_no)';

% Define exponential functions for one-step Richardson extrapolation 
x1 = 2^(-20);
Iz1 = exp(-x1*xip_z.^2);
Iz2 = exp(-(x1/2)*xip_z.^2);
Ir1 = exp(-x1*xip_r.^2);
Ir2 = exp(-(x1/2)*xip_r.^2);

% -------------------------------------------------------------------------
% DISPLACEMENTS
% -------------------------------------------------------------------------
% first part incl. J0: solved outside sum for all integration points 
% and then included inside sumation below

%%%%% u_z %%%%%
uzs = besselj(0,xip_z.*repmat(rho,1,nz)).*...
    (...
    (Aa-Cc.*(2-4*repmat(Nu,1,nz)-xip_z.*repmat(z/H,1,nz))).*exp(-xip_z.*repmat(Lami-z/H,1,nz))-...
    (Bb+Dd.*(2-4*repmat(Nu,1,nz)+xip_z.*repmat(z/H,1,nz))).*exp(-xip_z.*repmat(z/H-Lami1,1,nz))...
    ); 

% Improve the convergence for points residing close to the surface
Iuz1 = q.*alpha.*(-1)*H.*(1+Nu)./Ee.*sum(uzs.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z.*Iz1,2);
Iuz2 = q.*alpha.*(-1)*H.*(1+Nu)./Ee.*sum(uzs.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z.*Iz2,2);
uz = (4.*Iuz2-Iuz1)./3;

%%%%% u_r %%%%%
urs = besselj(1,xip_r.*repmat(rho,1,nz)).*...
    (...
    (Aa+Cc.*(1+xip_r.*repmat(z/H,1,nz))).*exp(-xip_r.*repmat(Lami-z/H,1,nz))+...
    (Bb-Dd.*(1-xip_r.*repmat(z/H,1,nz))).*exp(-xip_r.*repmat(z/H-Lami1,1,nz))...
    );

% Improve the convergence for points residing close to the surface
Iur1 = q.*alpha.*H.*(1+Nu)./Ee.*sum(urs.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r.*Ir1,2);
Iur2 = q.*alpha.*H.*(1+Nu)./Ee.*sum(urs.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r.*Ir2,2);
ur = (4*Iur2-Iur1)/3;

% -------------------------------------------------------------------------
% STRESSES 
% -------------------------------------------------------------------------
%%%%% Sigma_z %%%%%
sigma_zs = -xip_z.*besselj(0,xip_z.*repmat(rho,1,nz)).*...
    (...
    (Aa-Cc.*(1-2*repmat(Nu,1,nz)-xip_z.*repmat(z/H,1,nz))).*exp(-xip_z.*repmat(Lami-z/H,1,nz))+...
    (Bb+Dd.*(1-2*repmat(Nu,1,nz)+xip_z.*repmat(z/H,1,nz))).*exp(-xip_z.*repmat(z/H-Lami1,1,nz))...
    );

% Improve the convergence for points residing close to the surface
Isz1    = q.*alpha.*sum(sigma_zs.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z.*Iz1,2);
Isz2    = q.*alpha.*sum(sigma_zs.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z.*Iz2,2);
sigma_z = (4*Isz2-Isz1)/3;

%%%%% sigma_r %%%%%
sigma_rs = (xip_r.*besselj(0,xip_r.*repmat(rho,1,nz))...
    -((besselj(1,xip_r.*repmat(rho,1,nz)))./repmat(rho,1,nz))).*...
    (...
    (Aa+Cc.*(1+xip_r.*repmat(z/H,1,nz))).*exp(-xip_r.*repmat(Lami-z/H,1,nz))+...
    (Bb-Dd.*(1-xip_r.*repmat(z/H,1,nz))).*exp(-xip_r.*repmat(z/H-Lami1,1,nz))...
    )+ 2*repmat(Nu,1,nz).*xip_r.*besselj(0,xip_r.*repmat(rho,1,nz)).*...
    (Cc.*exp(-xip_r.*repmat(Lami-z/H,1,nz))-Dd.*exp(-xip_r.*repmat(z/H-Lami1,1,nz)));

% Improve the convergence for points residing close to the surface
Isr1    = q.*alpha.*sum(sigma_rs.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r.*Ir1,2);
Isr2    = q.*alpha.*sum(sigma_rs.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r.*Ir2,2);
sigma_r = (4*Isr2-Isr1)/3;

%%%%% sigma_theta %%%%%
sigma_thetas = (besselj(1,xip_r.*repmat(rho,1,nz))./repmat(rho,1,nz)).*...
    (...
    (Aa+Cc.*(1+xip_r.*repmat(z/H,1,nz))).*exp(-xip_r.*repmat(Lami-z/H,1,nz))+...
    (Bb-Dd.*(1-xip_r.*repmat(z/H,1,nz))).*exp(-xip_r.*repmat(z/H-Lami1,1,nz))...
    )+ 2*repmat(Nu,1,nz).*xip_r.*besselj(0,xip_r.*repmat(rho,1,nz)).*...
    (Cc.*exp(-xip_r.*repmat(Lami-z/H,1,nz))-Dd.*exp(-xip_r.*repmat(z/H-Lami1,1,nz)));

% Improve the convergence for points residing close to the surface
Isr1    = q.*alpha.*sum(sigma_thetas.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r.*Ir1,2);
Isr2    = q.*alpha.*sum(sigma_thetas.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r.*Ir2,2);
sigma_theta = (4*Isr2-Isr1)/3;

%%%%% tau_rz %%%%%
tau_rzs = xip_z.*besselj(1,xip_z.*repmat(rho,1,nz)).*...
    (...
    (Aa+Cc.*(2*repmat(Nu,1,nz)+xip_z.*repmat(z/H,1,nz))).*exp(-xip_z.*repmat(Lami-z/H,1,nz))-...
    (Bb-Dd.*(2*repmat(Nu,1,nz)-xip_z.*repmat(z/H,1,nz))).*exp(-xip_z.*repmat(z/H-Lami1,1,nz))...
    );

% Improve the convergence for points residing close to the surface
Itr1    = q.*alpha.*sum(tau_rzs.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z.*Iz1,2);
Itr2    = q.*alpha.*sum(tau_rzs.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z.*Iz2,2);
tau_rz = (4*Itr2-Itr1)/3;

% -------------------------------------------------------------------------
% Transformation between (r,theta) and (x,y) coordinates
% -------------------------------------------------------------------------
% Below a transformation of each response (evaluated above from
% the (r,theta) coordinates to the (x,y) coordinates are made. After the
% transformation, the responses are added together.

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
SIN = (Xd(:,2)-Xl(:,2))./r; % Minus on y-coordinate because z-axis points
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

% Organize stress matrix
sig1 = sparse(length(r)*3,length(r)*3,(3*3-4)*length(r));
sig1(1:3*3*length(r)+3:end)               = sigma_r;
sig1(3*length(r)+2:3*3*length(r)+3:end)   = sigma_theta;
sig1(2*3*length(r)+3:3*3*length(r)+3:end) = sigma_z;
sig1(3:3*3*length(r)+3:end)               = tau_rz;
sig1(2*3*length(r)+1:3*3*length(r)+3:end) = tau_rz;

% Transform deformations in to (x,y,z)-coordinates
u1 = T'*u1;
u  = u1;

% Transform stresses into (x,y,z)-coordinates
sig1 = T'*sig1*T;
sigm  = sig1;

% Transform 3x3 matrices into 6x1 vectors
dos            = zeros(6*xd*xl,1);
dos(1:6:end-5) = 1:9*xd*xl+3:9*xd*xl*length(r);
dos(2:6:end-4) = dos(1:6:end-5) + 3*xd*xl+1;
dos(3:6:end-3) = dos(1:6:end-5) + 2*3*xd*xl+2;
dos(4:6:end-2) = dos(1:6:end-5) + 1;
dos(5:6:end-1) = dos(2:6:end-4) + 1;
dos(6:6:end)   = dos(1:6:end-5) + 2;

sigv = sigm(dos);

% Organize matrix Txyz for adding displacements together
unos                = zeros(3*xd*xl,1);      % Number of ones in Txyz
unos(1:xl)          = 1:9*xd:9*xl*xd;
unos(xl+1:2*xl)     = unos(1:xl)+3*xd+1;
unos(2*xl+1:3*xl)   = unos(1:xl)+2*3*xd+2;

for i=1:xd-1
    unos(3*xl+1+3*xl*(i-1):3*xl+3*xl*i) = unos(1:3*xl)+i*(3*xd*3*xl+3);
end

% Matrix adding contributions together from different loads
Txyz       = spalloc(3*xd,3*length(r),length(unos)); 
Txyz(unos) = 1;

% Organize matrix Mxyz for adding displacements together
tres                = zeros(6*xd*xl,1);        % Number of ones in Txyz
tres(1:xl)          = 1:9*4*xd:9*4*xl*xd;
tres(xl+1:2*xl)     = tres(1:xl)+6*xd+1;
tres(2*xl+1:3*xl)   = tres(1:xl)+2*6*xd+2;
tres(3*xl+1:4*xl)   = tres(1:xl)+3*6*xd+3;
tres(4*xl+1:5*xl)   = tres(1:xl)+4*6*xd+4;
tres(5*xl+1:6*xl)   = tres(1:xl)+5*6*xd+5;

for i=1:xd-1
    tres(6*xl+1+6*xl*(i-1):6*xl+6*xl*i)=tres(1:6*xl)+i*(6*xd*6*xl+6);
end

% Matrix adding contributions together from different loads
Mxyz = spalloc(6*xd,6*length(r),length(tres));
Mxyz(tres) = 1;

% Final response vector

% Displacements
u          = Txyz*u;
alva.ux    = u(1:3:end-2);
alva.uy    = u(2:3:end-1);
alva.uz    = u(3:3:end);  

% Stresses
sig = -Mxyz*full(sigv);     % minus inserted to fullfil inwards positive 
                            % sign convention
alva.sigx  = sig(1:6:end-5); 
alva.sigy  = sig(2:6:end-4); 
alva.sigz  = sig(3:6:end-3); 
alva.sigxy = sig(4:6:end-2); 
alva.sigyz = sig(5:6:end-1); 
alva.sigxz = sig(6:6:end);   

% Strains
Eel = Ee(1:xl:end); % Reduce vector to evaluation points only
Nul = Nu(1:xl:end); % Reduce vector to evaluation points only

alva.epsx = 1./Eel.*(alva.sigx-Nul.*(alva.sigy+alva.sigz));
alva.epsy = 1./Eel.*(alva.sigy-Nul.*(alva.sigz+alva.sigx));
alva.epsz = 1./Eel.*(alva.sigz-Nul.*(alva.sigx+alva.sigy));
alva.epsxy = (1+Nul)./Eel.*alva.sigxy;
alva.epsyz = (1+Nul)./Eel.*alva.sigyz;
alva.epsxz = (1+Nul)./Eel.*alva.sigxz;
