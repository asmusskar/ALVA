function o = inv_loop(X,alva)

%--------------------------------------------------------------------------
% This file is an example file that calculates the displacements measured
% from a Falling Weight Deflectometer (FWD) test based on the layered elastic 
% theory and minimize error compared to measured values while predicting 
% E-modulus for each layer

% INPUT PARAMETERS
% E    : Layer Young's moduli
% Xd   : Sensor locations / evaluation points
%--------------------------------------------------------------------------

if length(X)==1
    E     = [X(1) X(1)];     % Initial Young's Modulus for single layer
else
    E     = X(1:1:end);
end

E      = (10.^(E.*8));

% Variable input parameters
alva.E = E;

% Constant input parameters
alva.Xd = [alva.gpos' zeros(length(alva.gpos),1) zeros(length(alva.gpos),1)];

% Calculate response
alva = init_LET(alva);

% Calculate (normalized) difference between measured and calculated 
% response
dd    = (alva.dFWD - alva.uz')./alva.dFWD;

% Error function (sum of squared relative differences)
o = 1/length(alva.dFWD)*sqrt(dd*dd'); 



