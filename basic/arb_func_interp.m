function ABCDi = arb_func_interp(n,xip,ABCD) 

%-------------------------------------------------------------------------
% DESCRIPTION:
% This function utilizes a cubic spline interpolation scheme to return 
% interpolated values for sought integration points (xip) based on the 
% calculated solutions for the coefficients of integration (arbitrary 
% functions).

% INPUT PARAMETERS:
% n:    Number of layers (including half space layer), i.e. n = length(E)
% xip:  Sought integration points
% ABCD: Arbitrary function values (soultion points) to be interpolated.
%-------------------------------------------------------------------------
m    = ABCD(end,:);
ABCD = ABCD(1:end-1,:);
for i=1:4*n
    ABCDi(i,:) = interp1(m,ABCD(i,:),xip,'PCHIP');
end