% DESCRIPTION:
% This script tests the implementation of the acclerated ALVA LET model, 
% calculating the surface displacements with length for a multilayered 
% pavement subjected to two circular loads utilizing the method proposed in 
% (Andersen, Levenberg, & Andersen, 2020). The results are compared to the 
% computer programme ELLEA1.
% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
% Andersen, S., Levenberg, E., & Andersen, M. B (2020). Efficient 
% reevaluation of surface displacements in a layered elastic half-space. 
% The International Journal of Pavement Engineering 21(4), 1-8. 
% https://doi.org/10.1080/10298436.2018.1483502
%
% Levenberg, E. (2016a). ELLEA1: Isotropic layered elasticity in excel: 
% Pavement analysis tool for students and engineers.

clear all, close all, clc

% -------------------------------------------------------------------------
% Select response analysis type
% -------------------------------------------------------------------------
alva.analysis = 'PolFit';
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
% Numerical parameters*
% -------------------------------------------------------------------------
alva.N  = 30;  % Number of Bessel zero points in numerical integration
alva.n  = 15;  % Number of Gauss points points between zero points.

% *Note: These should be small compared to conventional analysis to save
% computational time

% -------------------------------------------------------------------------
% Pavement material properties (minimum two layers required)
% -------------------------------------------------------------------------
alva.zi = [150 750];         % Depth of first n-1 layers from the surface 
                             % [mm]: last z = inf, and should not be added. 
                             % NB: zi(i) > zi(i-1) > z(i-2) ...
alva.E  = [3000 200 40];     % Layer Young's moduli [MPa]
alva.nu = [0.30 0.35 0.4];   % Layer Poisson's ratio [-]
alva.kh = [1e9 1e9];         % Interface bonding/horizontal spring [MPa/mm]

% -------------------------------------------------------------------------
% Load configuration
% -------------------------------------------------------------------------
alva.q  = [1.1
           1.1];         % Load pressure [MPa] (uniform vertical pressure)
alva.a  = [150
           150];         % Load radii [mm] (circular load)
alva.Xl = [0.0 -150
           0.0  150];    % Load positions [mm]: [x1 y1; x2 y2;..xi yi]; 

% -------------------------------------------------------------------------
% Location of evaluation points: [x1 y1 z1; x2 y2 z2;..]
% -------------------------------------------------------------------------
alva.Xd = [0    0   0;   50   0   0;  75   0   0;  100	0	0;  150	0	0;    
         200    0	0;  250   0   0; 300   0   0;  350	0	0;  400	0	0; 
         450    0	0;  500   0   0; 750   0   0; 1000	0	0; 1500	0	0;
        2000    0	0; 3000   0   0];
   
% -------------------------------------------------------------------------
% Initialize system and get response
% -------------------------------------------------------------------------

alva = init_LET(alva);

% Displacements
ux    = alva.ux;   % [mm]
uy    = alva.uy;   % [mm]
uz    = alva.uz;   % [mm]

% -------------------------------------------------------------------------
% Validation
% -------------------------------------------------------------------------

% Validation using independent software
alva.validation = 'duallength';         % select dual wheel load and 
                                        % displacements along x-axis at 
                                        % surface       
ellea           = validation_let(alva); % Collumns: [length/depth,sigx
                                        % ,sigy,sigx,sigzy,sigzx,sigxy,ux,
                                        % uy,uz]  
% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------

subplot(2,1,1), title('x-axis displacement (y=z=0)');
hold on, grid on
plot(alva.Xd(:,1),ux,':ks','LineWidth',1.25,'MarkerSize',5)
hold on
plot(ellea(:,1),ellea(:,8),':b+','LineWidth',1,'MarkerSize',5)
hold on
set(gca,'FontSize',9)
legend({'ALVA','ELLEA1'},'Location','SouthEast')
axis ij
xlabel('Distance from load, x-axis [mm]')
ylabel('Displacement, u_{x} [mm]')
hold on

subplot(2,1,2), title('z-axis displacement  (y=z=0)');
hold on, grid on
plot(alva.Xd(:,1),uz,':ks','LineWidth',1.25,'MarkerSize',5)
hold on
plot(ellea(:,1),ellea(:,10),':b+','LineWidth',1,'MarkerSize',5)
hold on
set(gca,'FontSize',9)
legend({'ALVA','ELLEA1'},'Location','SouthEast')
axis ij
xlabel('Distance from load, x-axis [mm]')
ylabel('Displacement, u_{z} [mm]')
hold off
