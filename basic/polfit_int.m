function I = polfit_int(Xd,Xl,a,H,Lam1,XipWip,Nu)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function evaluates the integrals Iz(i) and Ir(i) (i = 1,..,4)
% in the approximation method used (see [1]).

% INPUT PARAMETERS:
% Xd    : Coordinates of evaluation points (x,y,z)
% Xl    : Coordinates of load (x,y)
% a     : Load radii vector
% H     : Depth of layer n-1 out of n layers (i.e. thicness of structure)
% Lam1  : Relative depth of the top layer
% XipWip: Integration points and weights
% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
%[1] Andersen, S., Levenberg, E., & Andersen, M. B (2020). Efficient 
%    reevaluation of surface displacements in a layered elastic half-space. 
%    The International Journal of Pavement Engineering 21(4), 1-8. 
%    https://doi.org/10.1080/10298436.2018.1483502
%--------------------------------------------------------------------------

% Organize alpha and rho required for the integration

% Number of loads (xl) and deformation points (xd)
xl = size(Xl,1); % Number of load points
xd = size(Xd,1); % Number of deformation points

% Extend number of deformation points
Xl    = repmat(Xl,[xd,1]);  % Extends matrix to number of evaluation points
Xd    = repmat(Xd,[xl,1]);  % Xd is reorganized below to correspond to the 
                            % correct Xl content/sequence

 % Reorganize Xd so the first xl rows correspond to the same deformation
% point, and following xl rows correspond to the next deformation point
% etc... see description of the vector r nine lines below to understand the 
% organization
                            
% Built dx vector used to reorganize/sort Xd in the right order
dx = repmat(1:xd,xl,1);
dx = dx(:)' + repmat(0:xd:xd*(xl-1),1,xd);
Xd = Xd(dx,:);

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

r = sqrt((Xd(:,1)-Xl(:,1)).^2 + (Xd(:,2)-Xl(:,2)).^2);

% Introduce the parameter rho (notice: since r is a vector, rho is too)
rho = r/H;

% Alpha - repeated by the number of deformation points all loads influence
alpha = repmat(a/H,[xd,1]);
%------------------------------------------------------------------------
% Integration points and weights

lrho  = length(rho);                % Number of load-to-displacement-points
xip_z = XipWip(0*lrho+1:1*lrho,:);
wip_z = XipWip(1*lrho+1:2*lrho,:);
xip_r = XipWip(2*lrho+1:3*lrho,:);
wip_r = XipWip(3*lrho+1:4*lrho,:);

% Number of columns in reorganized integration points and weights (= number
% of integration points per integral)
nz = size(xip_z,2); 

% Perform precalculated integrals of polynomial approximations

% First integral Iz1
Iz1 = besselj(0,xip_z.*repmat(rho,1,nz)).*(xip_z.^2.*exp(-2*Lam1*xip_z));
     
Iz1 = alpha.*(-1)*H.*(1+Nu).*...
      sum(Iz1.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z,2);

% Second integral, Iz2
Iz2 = besselj(0,xip_z.*repmat(rho,1,nz)).*(xip_z.*exp(-2*Lam1*xip_z));
     
Iz2 = alpha.*(-1)*H.*(1+Nu).*...
      sum(Iz2.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z,2);

% Third integral, Iz3
Iz3 = besselj(0,xip_z.*repmat(rho,1,nz)).*exp(-2*Lam1*xip_z);
     
Iz3 = alpha.*(-1)*H.*(1+Nu).*...
      sum(Iz3.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z,2);
  
% Fourth integral, Iz4
Iz4 = besselj(0,xip_z.*repmat(rho,1,nz)).*(-(2-2*Nu));
     
Iz4 = alpha.*(-1)*H.*(1+Nu).*...
      sum(Iz4.*besselj(1,xip_z.*repmat(alpha,1,nz))./xip_z.*wip_z,2);

% First integral, Ir1
Ir1 = besselj(1,xip_r.*repmat(rho,1,nz)).*(xip_r.^2.*exp(-2*Lam1*xip_r));   

Ir1 = alpha.*H*(1+Nu).*sum(Ir1.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r,2);  

% Second integral, Ir2
Ir2 = besselj(1,xip_r.*repmat(rho,1,nz)).*(xip_r.*exp(-2*Lam1*xip_r));

Ir2 = alpha.*H*(1+Nu).*sum(Ir2.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r,2);  

% Third integral, Ir3
Ir3 = besselj(1,xip_r.*repmat(rho,1,nz)).*exp(-2*Lam1*xip_r);

Ir3 = alpha.*H*(1+Nu).*sum(Ir3.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r,2);  

% Foruth integral, Ir4
Ir4 = besselj(1,xip_r.*repmat(rho,1,nz)).*(-(1-2*Nu));

Ir4 = alpha.*H*(1+Nu).*sum(Ir4.*besselj(1,xip_r.*repmat(alpha,1,nz))./xip_r.*wip_r,2);  

% Save all integration values in one matrix
I = [Iz1 , Iz2 , Iz3 , Iz4 , Ir1 , Ir2 , Ir3 , Ir4];