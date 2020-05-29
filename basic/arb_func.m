function ABCD = arb_func(n,zi,E,nu,alva) 
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function evaluates the coefficients of integration Ai, Bi, Ci and Di
% of each layer. These are unitless functions that embody the layered 
% system properties and connectivity.

% INPUT PARAMETERS
% n    : Number of layers (including half space layer)
% zi   : Distance from surface to the bottom of each layer
% E    : Layer Young's moduli
% nu   : Layer Poissons ratio's
% bond : Interface bonding type
%--------------------------------------------------------------------------

% In order to speed the computational time, the number of matrix inversions
% is limited to 96, corresponding to 96 predetermined values of the
% integration variable m in the range of 0 to 100,000, as follows:
m = [1e-10, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60,...
    1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80,...
    4.00, 4.20, 4.40, 4.60, 4.80, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50,...
    8.00, 8.50, 9.00, 9.50, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00,...
    16.00, 17.00, 18.00, 19.00, 20.00, 25.00, 30.00, 35.00, 40.00,...
    45.00, 50.00, 55.00, 60.00, 65.00, 70.00, 75.00, 80.00, 85.00,...
    90.00, 95.00, 100.00, 110.00, 120.00, 130.00, 140.00, 150.00,...
    160.00, 170.00, 180.00, 190.00, 200.00, 210.00, 220.00, 230.00,...
    240.00, 250.00, 260.00, 270.00, 280.00, 290.00, 300.00, 350.00,...
    400.00, 450.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00,...
    2000.00, 10000.00, 100000.00];

% Number of integration points
Lm   = length(m);               % Lm = Length of m
ONES = ones(1,Lm); 

% Evaluate the parameters ABCD
H     = zi(end);                % Depth to the upper boundary of the 
                                % inifnite (bottom) layer from the top 
                                % surface 
lam   = [0 zi max(zi)*1e20]'/H; % Number of lambda values correspond to 
                                % the number of layer interfaces (including
                                % the surface of the top layer). We have 
                                % added a infinite value (1e20*max(zi)) in 
                                % the last place. But it will not influence 
                                % on the outcome, since it will only be 
                                % used in the last A and C coefficients 
                                % equations which are remmoved from the 
                                % system of equations before solved, since 
                                % these A and C values are zero. These two 
                                % equations are however organized, as it is 
                                % easier to do in the loop instead of 
                                % having a special case at the end.

% lam = [lam0 lam1 lam2 ... lam_n-1 0] = [z0 z1 z2 ... z_n-1 0]/H
% I.e. lami = lam(i+1)

% Make sure that E and nu are "standing" vectors in order for the
% multiplications below to be correct
if size(E,1) < size(E,2)
   E = E';
end

if size(nu,1) < size(nu,2)
   nu = nu';
end

if length(zi) > length(E)
    display('The length of zi is too high or the length of E is to short')
end
%--------------------------------------------------------------------------
% Setup the boundary and continuity conditions
%--------------------------------------------------------------------------
% Define the i and j-coordinates
j = 1:n-1;
i = j+1;

% Define the equations

% Equations based on the vertical stresses
sigmaz = [                                                    ... 
 1*ones(n-1,1)*ONES                                           ...  % Ai
 exp(-(lam(i)-lam(i-1))*m)                                    ...  % Bi
-(1-2*nu(j)*ONES-lam(i)*m)                                    ...  % Ci
 (1-2*nu(j)*ONES+lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)         ...  % Di
-exp(-(lam(i+1)-lam(i))*m)                                    ...  % A{i+1}
-1*ones(n-1,1)*ONES                                           ...  % B{i+1}
 (1-2*nu(j+1)*ONES-lam(i)*m).*exp(-(lam(i+1)-lam(i))*m)       ...  % C{i+1}
-(1-2*nu(j+1)*ONES+lam(i)*m)                                  ...  % D{i+1}
];

% Equations based on the shear stresses
taurz  = [                                                     ... 
 1*ones(n-1,1)*ONES                                            ... % Ai
-exp(-(lam(i)-lam(i-1))*m)                                     ... % Bi
 (2*nu(j)*ONES+lam(i)*m)                                       ... % Ci
 (2*nu(j)*ONES-lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)            ... % Di
-exp(-(lam(i+1)-lam(i))*m)                                     ... % A{i+1}
 1*ones(n-1,1)*ONES                                            ... % B{i+1}
-(2*nu(j+1)*ONES+lam(i)*m).*exp(-(lam(i+1)-lam(i))*m)          ... % C{i+1}
-(2*nu(j+1)*ONES-lam(i)*m)                                     ... % D{i+1}
];

if strcmp(alva.bond,'Frictionless')
    taurz1  = [                                                ...
     1*ones(n-1,1)*ONES                                        ... % Ai
    -exp(-(lam(i)-lam(i-1))*m)                                 ... % Bi
     (2*nu(j)*ONES+lam(i)*m)                                   ... % Ci
     (2*nu(j)*ONES-lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)        ... % Di
     -exp(-(lam(i+1)-lam(i))*m).*0                             ... % A{i+1}
      1*ones(n-1,1)*ONES.*0                                    ... % B{i+1}
     -(2*nu(j+1)*ONES+lam(i)*m).*exp(-(lam(i+1)-lam(i))*m).*0  ... % C{i+1}
     -(2*nu(j+1)*ONES-lam(i)*m).*0                             ... % D{i+1}
     ];
   taurz2  = [                                                 ...
     1*ones(n-1,1)*ONES.*0                                     ... % Ai
    -exp(-(lam(i)-lam(i-1))*m).*0                              ... % Bi
     (2*nu(j)*ONES+lam(i)*m).*0                                ... % Ci
     (2*nu(j)*ONES-lam(i)*m).*exp(-(lam(i)-lam(i-1))*m).*0     ... % Di
     -exp(-(lam(i+1)-lam(i))*m)                                ... % A{i+1}
      1*ones(n-1,1)*ONES                                       ... % B{i+1}
     -(2*nu(j+1)*ONES+lam(i)*m).*exp(-(lam(i+1)-lam(i))*m)     ... % C{i+1}
     -(2*nu(j+1)*ONES-lam(i)*m)                                ... % D{i+1}
     ];     
end

% Equations based on the vertical displacements             
uz = [                                                                              ...
 (1+nu(j))./E(j)*ONES                                                               ... % Ai
-(1+nu(j))./E(j)*ONES.*exp(-(lam(i)-lam(i-1))*m)                                    ... % Bi
-(1+nu(j))./E(j)*ONES.*(2-4*nu(j)*ONES-lam(i)*m)                                    ... % Ci
-(1+nu(j))./E(j)*ONES.*(2-4*nu(j)*ONES+lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)         ... % Di
-(1+nu(j+1))./E(j+1)*ONES.*exp(-(lam(i+1)-lam(i))*m)                                ... % A{i+1}
 (1+nu(j+1))./E(j+1)*ONES                                                           ... % B{i+1}
 (1+nu(j+1))./E(j+1)*ONES.*(2-4*nu(j+1)*ONES-lam(i)*m).*exp(-(lam(i+1)-lam(i))*m)   ... % C{i+1}
 (1+nu(j+1))./E(j+1)*ONES.*(2-4*nu(j+1)*ONES+lam(i)*m)                              ... % D{i+1} 
];                     
               
% Equations based on the radial displacements
ur = [                                                             ...
(1+nu(j))./E(j)*ONES                                               ... % Ai
(1+nu(j))./E(j)*ONES.*exp(-(lam(i)-lam(i-1))*m)                    ... % Bi
(1+nu(j))./E(j)*ONES.*(1+lam(i)*m)                                 ... % Ci
-(1+nu(j))./E(j)*ONES.*(1-lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)     ... % Di
-(1+nu(j+1))./E(j+1)*ONES.*exp(-(lam(i+1)-lam(i))*m)               ... % A{i+1}
-(1+nu(j+1))./E(j+1)*ONES                                          ... % B{i+1}
-(1+nu(j+1))./E(j+1)*ONES.*(1+lam(i)*m).*exp(-(lam(i+1)-lam(i))*m) ... % C{i+1}
(1+nu(j+1))./E(j+1)*ONES.*(1-lam(i)*m)                             ... % D{i+1}  
];

if strcmp(alva.bond,'Slip')
    kh = alva.kh;
    if size(kh,1) < size(kh,2)
        kh = kh';
    end
    ur = [                                                                                               ...
    -kh(j).*(H./m).*(1+nu(j))./E(j)-1*ones(n-1,1)*ONES                                                   ...  % Ai
     (1*ones(n-1,1)*ONES - kh(j).*(H./m).*(1+nu(j))./E(j)).*exp(-(lam(i)-lam(i-1))*m)                    ...  % Bi
    -kh(j).*(H./m).*(1+nu(j))./E(j).*(1+lam(i)*m)-(2*nu(j)*ONES+lam(i)*m)                                ...  % Ci
     (kh(j).*(H./m).*(1+nu(j))./E(j).*(1-lam(i)*m)-(2*nu(j)*ONES-lam(i)*m)).*exp(-(lam(i)-lam(i-1))*m)   ...  % Di
     kh(j).*(H./m).*(1+nu(j+1))./E(j+1).*exp(-(lam(i+1)-lam(i))*m)                                       ...  % A{i+1}
     kh(j).*(H./m).*(1+nu(j+1))./E(j+1)                                                                  ...  % B{i+1}
     kh(j).*(H./m).*(1+nu(j+1))./E(j+1).*(lam(i)*m+1).*exp(-(lam(i+1)-lam(i))*m)                         ...  % C{i+1}
     kh(j).*(H./m).*(1+nu(j+1))./E(j+1).*(lam(i)*m-1)                                                    ...  % D{i+1}
     ];
end

% Arrange the set of linear matrix system to solve for the parameters Ai, 
% Bi, Ci and Di
% Surface boundary conditions
BC0 = [exp(-m*(lam(2)-lam(1)))   1*ONES -(1-2*nu(1))*exp(-m*(lam(2)-lam(1)))  (1-2*nu(1))*ONES spalloc(1,4*Lm,0); % lam(2) = lam1, nu(1) = nu1
       exp(-m*(lam(2)-lam(1)))  -1*ONES    2*nu(1)*exp(-m*(lam(2)-lam(1)))      2*nu(1)*ONES   spalloc(1,4*Lm,0)];

% All conditions in the intermediate layers
BCs = [sigmaz
       taurz
        uz
        ur   ];

if strcmp(alva.bond,'Frictionless')
    BCs = [sigmaz
           uz
           taurz1
           taurz2];
end

% When inserting m in order to evaluate BC0 and BCs (=BC), the first 
% Lm = length(m) columns of BC0 will be exp(-m.*...) with different m 
% values. The following Lm columns will contain 1 and -1. We want to
% reorganize such that the first columns contains the equations
% in BC0 with m(1), the following columns correspond to the equations 
% in BC0 with the m(2) value etc., i.e., we want to go from
% BC0 = [exp(-m(1).*...)  exp(-m(2).*...)  exp(-m(3).*...) .. ]
% to
% BC0 = [exp(-m(1).*...)  1  -(1-2*nu(1))*exp(-m(1)...) .. ]

% The same reorganization we want to do with the BCs matrix
% Vectors used to change the the column sequence are denoted indx_c

% indx(:)' = [0 0 0 .. 1 1 1 .. 2 2 2 .. Lm-1 Lm-1 Lm-1]: 
% 8 = number of equations (=columns) in BCs for one m-value
indx_c = repmat(0:Lm-1,8,1); 
indx_c = repmat(1:Lm:(8-1)*Lm+1,1,Lm) + indx_c(:)';

% Also the rows of BCs should be reorganized in order for the lambda, nu 
% and E-values to be in the correct sequence. The index vector for this is
% denoted indx_r and is defined as
indx_r = [1:n-1 ; n:2*(n-1) ; 2*n-1:3*(n-1) ; 3*n-2:4*(n-1)];

% Use indx_c and indx_r to reorganice BCs and BC0 and order them in a
% united matrix BC. Only the rows of BCs are to be reorganized. BC0 is fine.
BC = [BC0(:,indx_c); BCs(indx_r(:),indx_c)];

% Now we have a matrix BC that contains all the coefficients for each
% integration point. BC has 4*n-2 rows (corresponding to the number of
% unknowns) and 4*n (the number of unknowns plus 2 - later we reduce this 
% to 4*n-2) columns per integration point, i.e. 4*n*Lm columns in
% total. Each (4*n-2) x 4*n - referred to as a submatrix, that only 
% represents a single integration point. We want to organize all submatrices 
% in a diagonal matrix of dimension Lm*(4*n-2) x 4*n*Lm, with each submatrix 
% decoupled from the others. 

% Below and indx vector that inserts the first submatrix in BC into the
% right positions in BCg is organized. The way to use indx is in the 
% following way: BCg(indx) = BC(:)

% indx for a single submatrix - step 1 of 5
indx = repmat(1:(4*n-2),8,1)'; % Indices on all cells in a single matrix. 
                               % We have (4*n-2)*8 coefficients in a single 
                               % sub matrix to fill in into the global matrix

% Additional content to indx is given (step 2 of 5). Two vectors a and c are 
% defined
a = repmat(0:7,4*n-2,1);
c = repmat(0:4:4*(n-2),4,1);
c = repmat([zeros(2,1) ; c(:)],8,1);

% Add vectors a and c to indx (step 3 of 5) - the content of indx will then
% represent the indeces for the first submatrix in the BCs matrix for a 
% single submatrix (the first one only)
indx = indx(:) + (a(:)+c)*(4*n-2)*Lm;

% Now indx is expanded to consider all Lm submatrices. A vector b is
% organized (step 4 of 5)
b = repmat([0 1:Lm-1],length(indx),1);

% Organize indx for all submatrices (step 5 of 5)
indx = repmat(indx,Lm,1) + b(:)*(4*n*Lm+1)*(4*n-2);

% Insert BC content into BCg in the right positions
% tic
% BCg(indx) = BC; % corresponds to BCg(indx) = BC(:)
% time_BCg=toc

% Indices to subindeces
[I,J] = ind2sub([Lm*(4*n-2) 4*n*Lm],indx);
BCg   = sparse(I,J,BC,Lm*(4*n-2),4*n*Lm);

% We now need to remove the the columns corresponding to the coefficients 
% An and Cn from each submatrix in BCg. The vector for this is denoted
% indx_cr and is defined as
indx_cr = [4*n-3:4*n:4*n*Lm-3; 4*n-1:4*n:4*n*Lm-1];

% Remove the columns
BCg(:,indx_cr) = [];

% Now we need to organize the external load vector. (use repmat for this)

% Define right-hand side of equations
Rhs = [1; spalloc(4*n-3,1,1)];  
Rhs = repmat(Rhs,Lm,1);     

% Evaluate A, B, C and D coefficients from the equation BCg*[A B C D ...]'
% = Rhs
ABCD0 = BCg\Rhs;                                                                         

% Redefine ABCD so that each column represents an integration point. 
% Furthermore zero values are introduced for A and C (which are needed to
% generalize the evaluations of sigma and u's in a simple way)
% ABCD = reorganize size 
ABCD    = zeros(4*n-2,Lm);      % Organize zero vector with two additinal 
                                % inputs compared to ABCD0
ABCD(:) = ABCD0;
ABCD    = [ABCD(1:end-2,:)      % [A1, B1, C1, D1, A2, B2 , .... , D{n-1}]
           spalloc(1,Lm,0)      % An = 0
           ABCD(end-1,:)        % Bn 
           spalloc(1,Lm,0)      % Cn = 0
           ABCD(end,:)    ];    % Dn

% Output A, B, C and D for all layers with m point as last row. To be used
% for interpolation at later stage.
ABCD    = full(ABCD);
ABCD    = [ABCD; m]; 