function [Az,Bz,Cz,Dz,Ar,Br,Cr,Dr]=arb_func_polfit(E,nu,zi,xip_z,xip_r,lr,alva)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function interpolates intermediate values for A, B, C and D by
% linear interpolation. This is done by taking one full sequence of A, B, C
% and D values for one 'load-to-deformation-point'-case and interpolate the
% remaning cases based on these values.

% INPUT PARAMETERS:
% xip_z , xip_r     : integration points
% zi , ri           : row numbers where xip_z and xip_r contain the highest 
%                     end value
% A1, B1, C1 and D1 : coefficients for uz evaluated in xip_z(zi,:)
% A2, B2, C2 and D2 : coefficients for ur evaluated in xip_r(ri,:)
% lr                : number of 'load-to-deformation-point' radii
%--------------------------------------------------------------------------

% Identify the largest xip values (i.e. largest interval) in xip_z and
% xip_r, respectively.
[xmax , xiz] = max(xip_z(:,end));
[xmax , xir] = max(xip_r(:,end));

% Define the rows with the largest interval
xipz = xip_z(xiz,:);
xipr = xip_r(xir,:);

% Evaluate A, B, C and D in the points xipz
ABCD = arb_func_plain(length(E),xipz,zi,E,nu,alva);

Az = ABCD(1,:); % Layer 1 A values
Bz = ABCD(2,:); %     ... B ...
Cz = ABCD(3,:); %     ... C ...
Dz = ABCD(4,:); %     ... D ...

% Evaluate A, B, C and D in the points xipr
ABCD = arb_func_plain(length(E),xipr,zi,E,nu,alva);
    
Ar = ABCD(1,:); % Layer 1 A values
Br = ABCD(2,:); %     ... B ...
Cr = ABCD(3,:); %     ... C ...
Dr = ABCD(4,:); %     ... D ...

% Number of integration points (and A, B, C and D values) per radius
nz = size(xip_z,2);

% Evaluate Az, Bz, Cz and Dz by interpolation from xip, A, B, C and D
xip_z = xip_z';
Az    = spline(xipz,Az,xip_z(:));
Bz    = spline(xipz,Bz,xip_z(:));
Cz    = spline(xipz,Cz,xip_z(:));
Dz    = spline(xipz,Dz,xip_z(:));

% Evaluate Ar, Br, Cr and Dr by interpolation from xip, A, B, C and D
xip_r = xip_r';
Ar    = spline(xipr,Ar,xip_r(:));
Br    = spline(xipr,Br,xip_r(:));
Cr    = spline(xipr,Cr,xip_r(:));
Dr    = spline(xipr,Dr,xip_r(:));

% Reorganize content of Az, Bz, Cz and Dz.
Az = reshape(Az,[nz,lr])';
Bz = reshape(Bz,[nz,lr])';
Cz = reshape(Cz,[nz,lr])';
Dz = reshape(Dz,[nz,lr])';

% Reorganize content of Ar, Br, Cr and Dr.
Ar = reshape(Ar,[nz,lr])';
Br = reshape(Br,[nz,lr])';
Cr = reshape(Cr,[nz,lr])';
Dr = reshape(Dr,[nz,lr])';