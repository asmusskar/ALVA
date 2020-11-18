% DESCRIPTION:
% This script tests the implementation of the ALVA LET model for inferring
% layer moduli, or so called "backcalculation" of layer moduli, based on 
% Falling Weight Deflectometer (FWD) measurements.

clear all, close all, clc

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
alva.zi = [161 421];        % Depth of first n-1 layers from the
                            % surface [mm]: last z = inf, and should not
                            % be added NB: zi(i) > zi(i-1) > z(i-2)...
alva.nu = [0.35 0.35 0.35]; % Layer Poisson's ratio [-]
alva.kh = [1e2 1e2];        % Interface bonding/horizontal spring [MPa/mm]

% -------------------------------------------------------------------------
% Load configuration
% -------------------------------------------------------------------------
alva.q  = 0.7;    % Load pressure [MPa] (uniform vertical pressure)
alva.a  = 150;    % Load radii [mm] (circular load)
alva.Xl = [0 0];  % Load positions [mm]: [x1 y1; x2 y2;..xi yi];

% -------------------------------------------------------------------------
% Sensor locations along x-axis in [mm]
% -------------------------------------------------------------------------

alva.gpos = [0 200 300 450 600 900 1200 1500 1800]; 

% -------------------------------------------------------------------------
% Response measurements (vertical surface displacements, uz) in [mm]
% -------------------------------------------------------------------------

alva.dFWD = [298.9 244.2 220.07 175.0 138.0 97.3 72.2 57.4 47.6].*1e-3; 

% -------------------------------------------------------------------------
% Initialize optimization problem
% -------------------------------------------------------------------------

% Stop criteria
tolfun  = 1e-4;             % Object function stop criteria
tolvar  = 1e-4;             % Step size value stop criteria
tolfval = 1e4;              % Maximum function evaluations
tolit   = 1e3;              % Maximum iterations 
options = optimset('PlotFcns',@optimplotfval,'MaxIter',tolit,...
    'MaxFunEvals',tolfval,'TolFun',tolfun,'TolX',tolvar);

% Initial parameters
alva.E = [200 200 200];     % Layer Young's moduli [MPa]
E0     = log10(alva.E)./8;  % Transform and scale parameters
x0     = E0;                % Variable input

% Run optimization algorithm
[xmin,fval,exitflag,output] = fminsearch(@(x)inv_loop(x,alva),x0,options);

fprintf('Iterations %d: ',output.iterations)
fprintf('Relative error %f\n',fval)

% -------------------------------------------------------------------------
% Get response for optimimal predicted E-moduli
% -------------------------------------------------------------------------
alva.Xd = [
      0      0      0;   50      0     0;  100      0      0;
    150      0      0;  200      0     0;  300      0      0;
    450      0      0;  600      0     0;  900      0      0;
    1200     0      0; 1500      0	   0; 1800      0      0];
alva.E = 10.^(xmin.*8);
alva   = init_LET(alva);

% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------

figure, title('Deflection curve');
hold on, grid on
plot(alva.gpos,alva.dFWD,'bo','LineWidth',1.25,'MarkerSize',5)
hold on
plot(alva.Xd(:,1),alva.uz,'-.r','LineWidth',1.0,'MarkerSize',6)
hold on
set(gca,'FontSize',9)
legend({'Measured','Predicted'},'Location','SouthEast')
axis ij
xlabel('Distance from load, x-axis [mm]')
ylabel('Displacement, u_{z} [mm]')
hold off

