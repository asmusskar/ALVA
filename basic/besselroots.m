function x = besselroots(o,N,k)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function looks up the zeros roots of Bessel function of order 0 and 
% 1 of the first kind. Up to 1000 points is considered.

% INPUT PARAMETERS
% o : Order of the bessel function
% N : Number of Bessel roots in integration
% k : kind: 1 or 2 
%     1) Bessel's differential equation that are finite at the origin
%     2) Bessel differential equation that have a singularity at the origin
%--------------------------------------------------------------------------

if o == 0 && k == 1
    load B0r_table 
    x = B0r(1:N); % Bessel 0 roots (B0r)
elseif o == 1 && k == 1
    load B1r_table 
    x = B1r(1:N); % Bessel 1 roots (B1r)
else
    disp('Not supported')
end
