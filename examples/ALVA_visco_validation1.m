% DESCRIPTION:
% This script tests the implementation of the ALVA VE model, calculating 
% the displacements for a single evaluation point on the surface of 
% a multilayered pavement considering a single circular load. The results 
% are compared to the computer programme ELLVA1, see details in 
% (Levenberg, 2016c).
% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
% Levenberg, E. (2016b). ELLVA1: Isotropic layered viscoelasticity in 
% excel (moving load): Advanced pavement analysis tool for students and 
% engineers.
% 
% Levenberg, E. (2016c). Viscoelastic pavement modeling with a 
% spreadsheet. Proceedings of the Eighth International Conference on 
% Maintenance and Rehabilitation of Pavements (mairepav8), 746–755. 
% doi:10.3850/978-981-11-0449-7-132-cd

clear all, close all, clc
% -------------------------------------------------------------------------
% Select response analysis type
% -------------------------------------------------------------------------
alva.analysis = 'Full';
% 1) 'Full'  : Conventional full integration (no approximations) applied 
%              for evaluation of displacements and stresses
% 2) 'PolFit': Use polynomial fit technique to reduce integral length
%              according to Andersen  et al. (2018) for evaluation of 
%              surface displacements

% -------------------------------------------------------------------------
% Select interface 
% -------------------------------------------------------------------------
alva.bond = 'Bonded';
% 1) 'Bonded'       : Full bonding between layers
% 2) 'Slip'         : Interface bonding factor
% 3) 'Frictionless' : No bonding between layers

% -------------------------------------------------------------------------
% Numerical parameters
% -------------------------------------------------------------------------
alva.N  = 300;  % Number of Bessel zero points in numerical integration
alva.n  = 30;   % Number of Gauss points points between zero points.

% -------------------------------------------------------------------------
% Pavement material properties (minimum two layers required)
% -------------------------------------------------------------------------
alva.zi = [150 750];         % Depth of first n-1 layers from the surface 
                             % [mm]: last z = inf, and should not be added. 
                             % NB: zi(i) > zi(i-1) > z(i-2) ...
alva.E  = [3000 200 40];     % Layer Young's moduli [MPa]
alva.nu = [0.30 0.35 0.4];   % Layer Poisson's ratio [-]
alva.kh = [1e9 1e9];         % Interface bonding/horizontal spring [MPa/mm]

% Viscoelastic properties of asphaltlayer E(1)
D0   = 2.5e-5;    % Instantaneous or glassy compliance, [1/MPa]
Dinf = 1.0e-2;    % Long time equilibrium or rubbery compliance [1/MPa]
nD   = 0.35;      % Slope of the creep curve in the transient region [-]
tauD = 1e3;       % Retardation time of the material response [s]

% -------------------------------------------------------------------------
% Load configuration
% -------------------------------------------------------------------------
alva.q   = 1.1;                 % Tire pressure vector [MPa]
alva.a   = 150;                 % Load radius vector [mm]
alva.Xl  = [0 300];             % Load coordinates in [mm]

% -------------------------------------------------------------------------
% Location of evaluation and mesh  points: [x1 y1 z1; x2 y2 z2;..]
% -------------------------------------------------------------------------
rep  = [0 0 0];                   % [x y z]-coordinate of evaluation point
nels = 200;                       % Number of elements
Vkh  = 60;                        % Vehicle speed [km/h]
V    = Vkh/3.6*1e3;               % Vehicle speed [mm/s]
x0   = 5000;                      % Start / end of mesh [mm]
dx0  = 2*x0/nels;                 % Mesh increment [mm]
dt   = dx0/V;  alva.dt = dt;      % Time increment [s]
tt   = 2*x0/V; alva.tt = tt;      % Total travel time [s]
xx   = (0:dx0:x0);                % Mesh x-direction
yy   = rep(2)*ones(length(xx),1); % Mesh y-direction
zz   = rep(3)*ones(length(xx),1); % Mesh z-direction
Xr   = [xx' yy zz]; alva.Xr = Xr; % Full mesh matrix

% Knots for cubic spline interpolation
kx   = [0 x0^(2/16) x0^(3/16) x0^(4/16) x0^(5/16) x0^(6/16) x0^(7/16)...
        x0^(8/16) x0^(9/16) x0^(10/16) x0^(11/16) x0^(12/16) x0^(13/16)...
        x0^(14/16) x0^(15/16) x0^(16/16)]; 
Xd   = [kx; ones(1,length(kx))*rep(2); ones(1,length(kx))*rep(3)]'; 
alva.Xd = Xd; 

% -------------------------------------------------------------------------
% Determine creep compliance curve and relaxation modulus
% -------------------------------------------------------------------------
alva.ti = 12;          % Number increments on compliance curve
alva    = VE_moduli(D0,Dinf,tauD,nD,alva);

% -------------------------------------------------------------------------
% Calculate linear elastic response at different times
% -------------------------------------------------------------------------
alva = VE_response(alva.E,alva);

% -------------------------------------------------------------------------
% Simulate moving load
% -------------------------------------------------------------------------

% Select output response: displacements (dxm, dym or dzm), stresses (sigxm, 
% sigym, sigzm, sigxym, sigyzm or sigxzm), strains (epsxm, epsym, epszm, 
% epsxym, epsyzm or epsxzm)

letres = alva.dzm;                   % output response dz                                           
veres  = VE_simulation(letres,alva); % run simulation

% -------------------------------------------------------------------------
% Validation 
% -------------------------------------------------------------------------

% Validation using independent software
ellva = validation_ve; % collumns: [position time uz]

% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------
figure, title('z-axis displacement vs. vehicle position (y=300 mm)');
plot((-x0:dx0:x0),veres,':k','LineWidth',1.25,'MarkerSize',6)
hold on, grid on
plot(ellva(:,1),ellva(:,3),'-.b','LineWidth',1.25,'MarkerSize',6)
hold on,
set(gca, 'FontSize',9)
legend({'ALVA','ELLVA VD'},...
    'Location','SouthEast')
axis ij
xlabel('Load position relative to evaluation point [mm]')
ylabel('Displacement, u_z [mm]')
hold off
