% DESCRIPTION:
% This script tests the implementation of the ALVA VE model, calculating 
% the displacements for a single evaluation point on the surface of 
% a multilayered pavement considering a single circular load. The results 
% are compared to the computer programme ELLVA1 [1], see details in [2].
% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
%[1] Levenberg, E. (2018). ELLVA_VD: Isotropic Layered Viscoelasticity in
%    Excel: Analysis tool for interpretation of deflections measured with a 
%    moving load 
%[2] Levenberg, E. (2016). Viscoelastic Pavement Modeling with a 
%    Spreadsheet. In Eighth International Conference on Maintenance and 
%    Rehabilitation of Pavements (pp. 746-755), Singapore

clear all, close all, clc
% -------------------------------------------------------------------------
% Definition of pavement system
% -------------------------------------------------------------------------

% y-> x ------------------- Pavement surface ---------------------------- %
% |                              layer 1      ^      ^       ^
% v z                                         | z1   |       |
%                                             v      |       |
% ---------------------------------------------------|-------|-------------
%                                layer 2             | z2    |
%                                                    v       |
%------------------------------------------------------------|-------------
%                                                            .
%                                                            .
% -------------------------------------------------------------------------
%                                layer n-1                   .
% -----------------------------------------------------------|-------------
%                                layer n                     | infinity
%                                                            v

% -------------------------------------------------------------------------
% Definition of coordinate system
% -------------------------------------------------------------------------
%                      o Xl (load point)
%  (z points downwards into the soil)
%  y -> x            ^  (deformation vector: points towards the load point)
%  |                /
%  v z             o Xd (deformation point)

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
kx   = [0 dx0 2*dx0 3*dx0 4*dx0 5*dx0 6*dx0 8*dx0 10*dx0 15*dx0 20*dx0...
       30*dx0 40*dx0 60*dx0 80*dx0 100*dx0]; 
Xd   = [kx; ones(1,length(kx))*rep(2); ones(1,length(kx))*rep(3)]'; 
alva.Xd = Xd; 

% -------------------------------------------------------------------------
% Determine creep compliance curve and relaxation modulus
% -------------------------------------------------------------------------
alva.ti = 12;          % Number increments on compliance curve
alva    = VE_moduli(D0,Dinf,tauD,nD,alva);

% -------------------------------------------------------------------------
% Calculate linear elastic response
% -------------------------------------------------------------------------
alva = VE_response(alva.E,alva);

% -------------------------------------------------------------------------
% Simulate moving load
% -------------------------------------------------------------------------

% Select output response: displacements (dx, dy or dz), stresses (sigx, 
% sigy, sigz, sigxy, sigyz or sigxz), strains (epsx, epsy, epsz, epsxy, 
% epsyz or epsxz)

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
legend({'ALVA','ELLVA VD'},...
    'Location','SouthEast', 'FontSize',9)
axis ij
xlabel('Load position relative to evaluation point [mm]')
ylabel('Displacement, u_z [mm]')
hold off
