% DESCRIPTION:
% This script tests the implementation of the ALVA LET model calculating 
% vertical stresses and displacements with depth for a half-space subjected
% to a single circular load.
% The results are compared to the analytical Boussinesq solution and the 
% computer programme ELLEA1 [1].
% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
%[1] Levenberg, E. (2016). ELLEA1: Isotropic Layered Elasticity in Excel:
%    Pavement analysis tool for students and engineers

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
% (x,y)-plane:             Pavement surface     
% z-direction:             Pavement depth              
% Xl_i = (xl_i,yl_i):      Coordinates of load i on pavement surface
% Xd_i = (xd_i,yd_i,zd_i): Coordinates of deformation point i
%
%                | Load 1 (Xl_1)              |  Load 2 (Xl_2)
%                v                            v
% y-> x  -----^-----------Pavement surface --------------------------------
% |            \                             ^
% v z    ^      o  Xd_1 (deformation point 1) \
%       /                                      o Xd_3 (deformation point 3)
%      0 Xd_2 (deformation point 2)

% -------------------------------------------------------------------------
% Select response analysis type
% -------------------------------------------------------------------------
alva.analysis = 'Full';
% 1) 'Full'     : Conventional full integration with one-step Richardson
%                 extrapolation (for improved convergence near surface) 
%                 applied for evaluation of displacements and stresses
% 2) 'PolFit'   : Use polynomial fit technique to reduce integral length
%                 according to Andersen et al. (2018) for evaluation of
%                 surface displacements

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
alva.zi = [150 750];         % Depth of first n-1 layers from the 
                             % surface [mm]: last z = inf, and should not 
                             % be added NB: zi(i) > zi(i-1) > z(i-2)...
alva.E  = [200 200 200];     % Layer Young's moduli [MPa]
alva.nu = [0.35 0.35 0.35];  % Layer Poisson's ratio [-]
alva.kh = [1e9 1e9];         % Interface bonding/horizontal spring [MPa/mm]

% -------------------------------------------------------------------------
% Load configuration
% -------------------------------------------------------------------------
alva.q  = 1.1;    % Load pressure [MPa] (uniform vertical pressure)
alva.a  = 150;    % Load radii [mm] (circular load)
alva.Xl = [0 0];  % Load positions [m]: [x1 y1; x2 y2;..xi yi];

% -------------------------------------------------------------------------
% Location of evaluation points: [x1 y1 z1; x2 y2 z2;..] in [mm]
% -------------------------------------------------------------------------
alva.Xd = [0      0      0; 0      0     10; 0      0      20;    
           0      0     30; 0      0     40; 0      0      50;	
           0      0     75; 0      0	100; 0      0     150;
           0      0    200; 0      0	250; 0      0     300; 
           0      0    350; 0      0    400; 0      0     450;
           0      0    500; 0      0    750; 0      0    1000];
% -------------------------------------------------------------------------
% Initialize system and get response
% -------------------------------------------------------------------------

alva = init_LET(alva);

% Displacements
ux    = alva.ux;
uy    = alva.uy;   % [mm]
uz    = alva.uz;   % [mm]

% Stresses
sigx  = alva.sigx;  % [MPa]
sigy  = alva.sigy;  % [MPa]
sigz  = alva.sigz;  % [MPa]
sigxy = alva.sigxy; % [MPa]
sigyz = alva.sigyz; % [MPa]
sigxz = alva.sigxz; % [MPa]

% Strains
epsx  = alva.epsx.*1e6;  % [micro strain]
epsy  = alva.epsy.*1e6;  % [micro strain]
epsz  = alva.epsz.*1e6;  % [micro strain]
epsxy = alva.epsxy.*1e6; % [micro strain]
epsyz = alva.epsyz.*1e6; % [micro strain]
epsxz = alva.epsxz.*1e6; % [micro strain]

% -------------------------------------------------------------------------
% Validation
% -------------------------------------------------------------------------
% Validation using analytical solution
boussinesq      = validation_boussinesq(alva); % Collumns: [sigz,uz]

% Validation using independent software
alva.validation = 'halfspace';          % select half-space solution         
ellea           = validation_let(alva); % Collumns: [length/depth,sigx
                                        % ,sigy,sigx,sigzy,sigzx,sigxy,ux,
                                        % uy,uz] 
% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------

subplot(1,2,1), title('z-axis displacement (x=y=0)');
hold on, grid on
plot(uz,alva.Xd(:,3),':ks','LineWidth',1.25,'MarkerSize',5)
hold on, grid on
plot(ellea(:,10),ellea(:,1),':b+','LineWidth',1,'MarkerSize',5)
hold on
plot(boussinesq(:,2),alva.Xd(:,3),'-.r','LineWidth',1,'MarkerSize',5)
hold on
legend({'ALVA','ELLEA1','Analytical'},...
    'Location','SouthEast', 'FontSize',9)
axis ij
xlabel('Displacement, u_{z} [mm]')
ylabel('Depth, z-axis [mm]')
hold on

subplot(1,2,2), title('z-axis stress (x=y=0)');
hold on, grid on
plot(sigz,alva.Xd(:,3),':ks','LineWidth',1.25,'MarkerSize',5)
hold on
plot(ellea(:,4),ellea(:,1),':b+','LineWidth',1,'MarkerSize',5)
hold on
plot(boussinesq(:,1),alva.Xd(:,3),'-.r','LineWidth',1,'MarkerSize',5)
hold on
legend({'ALVA','ELLEA1','Analytical'},...
    'Location','SouthEast', 'FontSize',9)
axis ij
xlabel('Stress, \sigma_{z} [MPa]')
ylabel('Depth, z-axis [mm]')
hold off