function alva = VE_response(E,alva)
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function calculates the elastic response at selected locations
% prior to spline interpolation

% INPUT PARAMETERS
% x    : Coordinates of knots
% ti   : Time increments
% E    : Layer Young's moduli [MPa]
% Et   : Relaxation modulus [MPa]
% alva : System parameters (e.g., load config etc.)
%--------------------------------------------------------------------------

Xd = alva.Xd;
ti = alva.ti;
Et = alva.Et;

% Initialize response parameters
dx   = zeros(length(Xd),length(ti));  dy = dx; dz    = dx;
sigx = dx; sigy  = dx; sigz  = dx; sigxy = dx; sigyz = dx; sigxz = dx;
epsx = dx; epsy  = dx; epsz  = dx; epsxy = dx; epsyz = dx; epsxz = dx;

for i=1:ti
    E(1)     = Et(i);
    alva.E   = E;
    alva     = init_LET(alva);
    
    % Displacements
    dx(:,i) = alva.ux;        % [mm]
    dy(:,i) = alva.uy;        % [mm]
    dz(:,i) = alva.uz;        % [mm]
    
    if strcmp(alva.analysis,'Full')
        % Stresses
        sigx(:,i)  = alva.sigx;   % [MPa]
        sigy(:,i)  = alva.sigy;   % [MPa]
        sigz(:,i)  = alva.sigz;   % [MPa]
        sigxy(:,i) = alva.sigxy;  % [MPa]
        sigyz(:,i) = alva.sigyz;  % [MPa]
        sigxz(:,i) = alva.sigxz;  % [MPa]
        
        % Strains
        epsx(:,i)  = alva.epsx.*1e6;   % [micro strain]
        epsy(:,i)  = alva.epsy.*1e6;   % [micro strain]
        epsz(:,i)  = alva.epsz.*1e6;   % [micro strain]
        epsxy(:,i) = alva.epsxy.*1e6;  % [micro strain]
        epsyz(:,i) = alva.epsyz.*1e6;  % [micro strain]
        epsxz(:,i) = alva.epsxz.*1e6;  % [micro strain]
    end
end

% Store displacements
alva.dxm = dx; alva.dym = dy; alva.dzm = dz;

% Store stresses
alva.sigxm  = sigx;  alva.sigym  = sigy;  alva.sigzm  = sigz;
alva.sigxym = sigxy; alva.sigyzm = sigyz; alva.sigxzm = sigxz;

% Store strains
alva.epsxm  = epsx;  alva.epsym  = epsy;  alva.epszm  = epsz;
alva.epsxym = epsxy; alva.epsyzm = epsyz; alva.epsxzm = epsxz;
