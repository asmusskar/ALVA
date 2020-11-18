% DESCRIPTION:
% This script tests the implementation of the ALVA LET model, calculating 
% responses at key locations for a multilayered pavement subjected to a 
% single load and a dual load respectively utilizing the method proposed. 
% The results are compared to benchmark results reported by the European 
% Commission (Belgian Road Research Centre, 2000).
% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
% Belgian Road Research Centre. (2000). Advanced models for analytical 
% design of european pavement structures. Brussels, Belgium: European 
% Commission. Retrieved from https://trimis.ec.europa.eu/project/...
% advanced-models-analytical-design-european-pavement-structures

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
% Numerical parameters
% -------------------------------------------------------------------------
alva.N  = 300;  % Number of Bessel zero points in numerical integration
alva.n  = 30;   % Number of Gauss points points between zero points.

% -------------------------------------------------------------------------
% Pavement material properties (minimum two layers required)
% -------------------------------------------------------------------------
alva.zi = [260 760];         % Depth of first n-1 layers from the 
                             % surface [mm]: last z = inf, and should not 
                             % be added NB: zi(i) > zi(i-1) > z(i-2)...
alva.E  = [5000 200 50];     % Layer Young's moduli [MPa]
alva.nu = [0.35 0.40 0.45];  % Layer Poisson's ratio [-]

% -------------------------------------------------------------------------
% Load configuration - Single case
% -------------------------------------------------------------------------

alva.q  = [0.7];         % Load pressure [MPa] (uniform vertical pressure)
alva.a  = [150.8];       % Load radii [mm] (circular load)
alva.Xl = [0.0 0.0];     % Load positions [mm]: [x1 y1; x2 y2;..xi yi];

% -------------------------------------------------------------------------
% Location of evaluation points: [x1 y1 z1; x2 y2 z2;..] - Single case
% -------------------------------------------------------------------------
alva.Xd = [    0 0 0;     0 0 259.9;     0 0 260.01;     0 0 760.1;   
           150.8 0 0; 151 0 259.9; 150.8 0 260.01; 150.8 0 760.1
           ];
       
% -------------------------------------------------------------------------
% Initialize system and get response - Single case
% -------------------------------------------------------------------------

% Case 1: Full bond both interfaces 
alva.bond = 'Bonded';
alva = init_LET(alva);

% Stresses
sigz1s  = alva.sigz;  % [MPa]

% Strains
epsx1s  = alva.epsx.*1e6;  % [micro strain]
epsy1s  = alva.epsy.*1e6;  % [micro strain]
epsz1s  = alva.epsz.*1e6;  % [micro strain]

% Case 2: Full bond both interfaces - kh is high 
alva.bond = 'Slip';
alva.kh = [1e6 1e6];      % Interface bonding/horizontal spring [MPa/mm]

alva = init_LET(alva);

% Stresses
sigz2s  = alva.sigz;  % [MPa]

% Strains
epsx2s  = alva.epsx.*1e6;  % [micro strain]
epsy2s  = alva.epsy.*1e6;  % [micro strain]
epsz2s  = alva.epsz.*1e6;  % [micro strain]


% Case 3: Reduceed bond interface 1 - kh is low
alva.bond = 'Slip';
alva.kh = [1e-6 1e6];         % Interface bonding/horizontal spring [MPa/mm]

alva = init_LET(alva);

% Stresses
sigz3s  = alva.sigz;  % [MPa]

% Strains
epsx3s  = alva.epsx.*1e6;  % [micro strain]
epsy3s  = alva.epsy.*1e6;  % [micro strain]
epsz3s  = alva.epsz.*1e6;  % [micro strain]

% -------------------------------------------------------------------------
% Load configuration - Dual case
% -------------------------------------------------------------------------

alva.q  = [0.7
           0.7];         % Load pressure [MPa] (uniform vertical pressure)
alva.a  = [106.6
           106.6];       % Load radii [mm] (circular load)
alva.Xl = [-170 0.0
            170 0.0];    % Load positions [mm]: [x1 y1; x2 y2;..xi yi];
          
% -------------------------------------------------------------------------
% Location of evaluation points: [x1 y1 z1; x2 y2 z2;..] - Dual case
% -------------------------------------------------------------------------

alva.Xd = [
           170 0 0; 170 0 259.99; 170 0 260.1; 170 0 760.1;
             0 0 0;   0 0 259.99;   0 0 260.1;   0 0 760.1
           ];   
       
% -------------------------------------------------------------------------
% Initialize system and get response
% -------------------------------------------------------------------------

% Case 1: Full bond both interfaces 
alva.bond = 'Bonded';
alva = init_LET(alva);

% Stresses
sigz1d  = alva.sigz;  % [MPa]

% Strains
epsx1d  = alva.epsx.*1e6;  % [micro strain]
epsy1d  = alva.epsy.*1e6;  % [micro strain]
epsz1d  = alva.epsz.*1e6;  % [micro strain]

% Case 2: Full bond both interfaces - kh is high 
alva.bond = 'Slip';
alva.kh = [1e6 1e6];      % Interface bonding/horizontal spring [MPa/mm]

alva = init_LET(alva);

% Stresses
sigz2d  = alva.sigz;  % [MPa]

% Strains
epsx2d  = alva.epsx.*1e6;  % [micro strain]
epsy2d  = alva.epsy.*1e6;  % [micro strain]
epsz2d  = alva.epsz.*1e6;  % [micro strain]

% Case 3: Reduceed bond interface 1 - kh is low
alva.bond = 'Slip';
alva.kh = [1e-6 1e6];         % Interface bonding/horizontal spring [MPa/mm]

alva = init_LET(alva);

% Stresses
sigz3d  = alva.sigz;  % [MPa]

% Strains
epsx3d  = alva.epsx.*1e6;  % [micro strain]
epsy3d  = alva.epsy.*1e6;  % [micro strain]
epsz3d  = alva.epsz.*1e6;  % [micro strain]

% -------------------------------------------------------------------------
% Validation: [sigmaz(c) epsh(c) epsz(c) epsz(c) sigmaz(e) epsh(e)
% epsz(e) epsz(e)]; where (c) and (e) indicate center and edge of load,
% respectively
% -------------------------------------------------------------------------

bond_single = [0.700 -100.5 251.7 185.0 0.350 -61.9 192.2 177.5;
               0.817 -100.5 251.6 185.3 0.319 -62.0 192.2 177.0;
               0.700 -100.5 251.6 185.1 0.345 -61.9 192.2 177.5];
           
slip_single = [0.700 -120  1 217  0.350  -78 -10 205;
               0.700 -120  1 216  0.000  -78 -10 205;
               0.700 -119 11 217  0.000  -77  -1 205];

bond_dual = [0.700     0   186  170   0.0000     0   182   177;
             1.466   -85   186  170  -0.0045   -89   183   177;
             0.700   -85   186  170   0.0000   -89   183   177];

slip_dual = [0.700    0   9  193 0.0000      0  -12  204;
             0.700 -120  -1  216 0.0000    -78  -10  205;
             0.700 -101  -3  194 0.0000   -106    1  205];


% -------------------------------------------------------------------------
% Plot results (table function not compatible with Octave)
% -------------------------------------------------------------------------

Location = {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8'};
Description     = [{'Vertical stress surface at center of load'};...
            {'Horizontal strain bottom layer 1 at center of load'};...
            {'Vertical strain top layer 2 at center of load'};...
            {'Vertical strain top layer 3 at center of load'};...
            {'Vertical stress surface at edge of load'};...
            {'Horizontal strain bottom layer 1 at edge of load'};...
            {'Vertical strain top layer 2 at edge of load'};...
            {'Vertical strain top layer 3 at edge of load'}];
%TableGuide = table(Description,'RowNames',Location)

TableGuide = [{'R1: Vertical stress surface at center of load'};...
            {'R2: Horizontal strain bottom layer 1 at center of load'};...
            {'R3: Vertical strain top layer 2 at center of load'};...
            {'R4: Vertical strain top layer 3 at center of load'};...
            {'R5: Vertical stress surface at edge of load'};...
            {'R6: Horizontal strain bottom layer 1 at edge of load'};...
            {'R7: Vertical strain top layer 2 at edge of load'};...
            {'R8: Vertical strain top layer 3 at edge of load'}]

CodeName1 = {'BISAR';'KENLAYER';'GAMES';'ALVA (slip)';'ALVA (bond)'}
CodeName2 = {'BISAR';'KENLAYER';'GAMES';'ALVA (slip)'};

% Single wheel bonded interfaces
A1sb = [bond_single(:,1);sigz2s(1);sigz1s(1)];
A2sb = [bond_single(:,2);epsx2s(2);epsx1s(2)];
A3sb = [bond_single(:,3);epsz2s(3);epsz1s(3)];
A4sb = [bond_single(:,4);epsz2s(4);epsz1s(4)];
A5sb = [bond_single(:,5);sigz2s(5);sigz1s(5)];
A6sb = [bond_single(:,6);epsx2s(6);epsx1s(6)];
A7sb = [bond_single(:,7);epsz2s(7);epsz1s(7)];
A8sb = [bond_single(:,8);epsz2s(8);epsz1s(8)];
%SingleWheelBond = table(A1sb,A2sb,A3sb,A4sb,A5sb,A6sb,A7sb,A8sb,...
%    'RowNames',CodeName1)
SingleWheelBond = [A1sb,A2sb,A3sb,A4sb,A5sb,A6sb,A7sb,A8sb]

% Single wheel unbonded interface 1, bonded interface 2
A1ss = [slip_single(:,1);sigz3s(1)];
A2ss = [slip_single(:,2);epsx3s(2)];
A3ss = [slip_single(:,3);epsz3s(3)];
A4ss = [slip_single(:,4);epsz3s(4)];
A5ss = [slip_single(:,5);sigz3s(5)];
A6ss = [slip_single(:,6);epsx3s(6)];
A7ss = [slip_single(:,7);epsz3s(7)];
A8ss = [slip_single(:,8);epsz3s(8)];
%SingleWheelSlip = table(A1ss,A2ss,A3ss,A4ss,A5ss,A6ss,A7ss,A8ss,...
%    'RowNames',CodeName2)
SingleWheelSlip = [A1ss,A2ss,A3ss,A4ss,A5ss,A6ss,A7ss,A8ss]

% Dual wheel bonded interfaces
A1db = [bond_dual(:,1);sigz2d(1);sigz1d(1)];
A2db = [bond_dual(:,2);epsy2d(2);epsy1d(2)];
A3db = [bond_dual(:,3);epsz2d(3);epsz1d(3)];
A4db = [bond_dual(:,4);epsz2d(4);epsz1d(4)];
A5db = [bond_dual(:,5);sigz2d(5);sigz1d(5)];
A6db = [bond_dual(:,6);epsy2d(6);epsy1d(6)];
A7db = [bond_dual(:,7);epsz2d(7);epsz1d(7)];
A8db = [bond_dual(:,8);epsz2d(8);epsz1d(8)];
%DualWheelBond = table(A1db,A2db,A3db,A4db,A5db,A6db,A7db,A8db,...
%    'RowNames',CodeName1)
DualWheelBond = [A1db,A2db,A3db,A4db,A5db,A6db,A7db,A8db]

% Dual wheel unbonded interface 1, bonded interface 2
A1ds = [slip_dual(:,1);sigz3d(1)];
A2ds = [slip_dual(:,2);epsy3d(2)];
A3ds = [slip_dual(:,3);epsz3d(3)];
A4ds = [slip_dual(:,4);epsz3d(4)];
A5ds = [slip_dual(:,5);sigz3d(5)];
A6ds = [slip_dual(:,6);epsy3d(6)];
A7ds = [slip_dual(:,7);epsz3d(7)];
A8ds = [slip_dual(:,8);epsz3d(8)];
%DualWheelSlip = table(A1ds,A2ds,A3ds,A4ds,A5ds,A6ds,A7ds,A8ds,...
%    'RowNames',CodeName2)
DualWheelSlip = [A1ds,A2ds,A3ds,A4ds,A5ds,A6ds,A7ds,A8ds]
