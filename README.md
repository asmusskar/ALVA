Skar, Asmus, Andersen, Sebastian and Julius Nielsen (2020): Adaptive Layered Viscoelastic Analysis (ALVA). Technical University of Denmark. Software. https://doi.org/10.11583/DTU.12387305

# ALVA
Adaptive Layered Viscoelastic Analysis (ALVA) is a MATLAB-package for pavement modeling. The aim with the software is to equip the civil engineering community with an advanced pavement modeling tool and computer package that is highly adaptive, transparent and open-access, capable of supporting current and future pavement evaluation needs. 

<div>
<img src="images/Roadmap_ALVA.png" width="90%">
</div>

The components of the ALVA software (shown in the Figure) is briefly described below:
* `main.m` - main script for defining vehicle loading conditions, pavement structure geometry and material properties, numerical parameters (for balancing accuracy and efficiency), analysis type and evaluation points (i.e., location of the output response), and post-processing of the results. `main.m` scripts can be found in the ../ALVA/examples_ folder.
* `init_LET.m` - initialize the analysis of a layered elastic half space selecting between: (a) `let_response_full.m` - evaluates the response (displacements, stresses and strains) of a layered elastic half-space model or (b) `let_response_polfit.m` - evaluates the displacements of a layered elastic half space model at the surface only. 
* `init_LET.m` also works as an interface application for simulating moving loads in the *Viscoelastic module*; based on the asphalt material properties from `VE_moduli.m` the module is engaged to calculate the elastic responses in `VE_response.m`, which is subsequently used to predict the linear viscoelastic response of a moving vehicle in `VE_simulation.m`.
* The main support functions in `init_LET.m` are (a) `arb_funct.m` - evaluates the coefficients of integration that embody the layered system properties and connectivity and (b) `numint_coeff.m` – organizes the integration points and weights for numerical integration for the user defined number of Bessel roots from `besselroots.m` and gauss points and weights between Bessel roots from `lookup_gauss.m`.

* Additional support functions for efficient surface displacement calculations in `let_response_polfit.m` is `arb_func_polfit.m` – interpolates intermediate coefficients of integration, `polfit_int.m` – evaluates coefficient proportional integrals, and `polfit_abc.m` – evaluates coefficients a, b and c for the polynomial.
Examples on ‘main.m’ scripts can be found in the _../ALVA/examples_ folder and all underlying components (i.e., source files) described above in the _../ALVA/basic_ folder.

## Main input parameters
The core algorithm behind this package is based on LET, i.e., the classic formulation for an <i>N</i>-layered half-space, shown in <b>Figure 1</b> 

<div>
<img src="images/N_layer.png" width="90%">
</div>
<p>
 <b>Figure 1:</b><i>N</i>-layered half-space model 
</p>

In this model all layers are assumed linear elastic, isotropic, homogeneous, fully bonded, and weightless. The model inputs include Young’s modulus E<sub>n</sub>, Poisson’s ratio <i>&Nu;<sub>n</sub></i>, and layer thickness t<sub>n</sub> (where <i>n</i> denotes the layer number). This model is engaged to calculate the response at any point, A<sub>j</sub>, of interest and for a given set of uniformly distributed circular loadings with load radius, <i>a</i>, and pressure <i>q</i>). An overview of LET model assumptions and solution procedure is given in [Khazanovich and Wang (2007)](https://journals.sagepub.com/doi/abs/10.3141/2037-06)

The viscoelastic response is approximated based on the LET calculations utilizing the methodology and load scheme suggested by [Levenberg (2016)](https://orbit.dtu.dk/en/publications/viscoelastic-pavement-modeling-with-a-spreadsheet) (see <b>Figure 2</b>). Viscoelastic layers are associated with a creep compliance and model parameters <i>D<sub>0</sub></i> and <i>D<sub>&infin;</sub></i>, the short and long time compliances (respectively), and shape parameters <i>&tau;<sub>D</sub></i> and <i>n<sub>D</sub></i>, controlling the transition between <i>D<sub>0</sub></i> and <i>D<sub>&infin;</sub></i>. 

<div>
<img src="images/VE_mesh.png" width="90%">
</div>
<p>
<b>Figure 2:</b>Load scheme to simulate movement
</p>

The load moves in a straight line from <i>x=-x<sub>0</sub></i> (Start) to <i>x=x<sub>0</sub></i> (End). The travel path is decomposed into <i>N</i> intervals (<i>i=1,…,N</i>), each <i>&Delta;x</i> long. The point of response evaluation <i>A<sub>j</sub></i> is indicated in the Figure; this point is located near the middle of the travel path (i.e., <i>x</i>-coordinate of zero), at <i>y</i>-coordinate <i>y<sub>0</sub></i> and depth <i>z<sub>0</sub></i> below the surface. 

## Validation examples
ALVA comes with six validation examples (i.e., `main.m` scripts) comparing ALVA to existing codes and analytical formulations. 
The example files can be found in the in the _../ALVA/examples_ folder and the results obtained with independent codes can be found in the found in the _../ALVA/validation_ folder.  
The following examples are included in this package:

* 'ALVA_bonding_validation1.m' - tests the implementation of the interface spring model, as well as compares the ALVA model with a range of commonly used softaware at critical positions within the pavament system published by the [European Commission](https://trimis.ec.europa.eu/project/advanced-models-analytical-design-european-pavement-structures).


| Layer    | Thickness (mm) | Youngs moduli (MPa) | Poisson's ratio  |    
|----------|----------------|---------------------|------------------|
| 1        |   260          |     5000            |   0.35           |
| 2        |   500          |      200            |   0.40           |
| 3        |   &infin;      |       50            |   0.45           |

<p>
<b>Table 1</b>: Reference pavement system used in basic validation.
</p>

| Load     | Radius (mm)    | Pressure, q (MPa) | x-position (mm) |  y-position (mm) |    
|----------|----------------|-------------------|-----------------|------------------|
| Single   |   150.8        |     0.7           |      0          |   0              |
| Dual (1) |   106.6        |     0.7           |   -170          |   0              |
| Dual (2) |   106.6        |     0.7           |    170          |   0              |
<p>
<b>Table 2</b>: Load cases used in basic validation. 
</p>

| Name        | Response                                   | Location                   | Unit   |
|-------------|--------------------------------------------|----------------------------|--------|
| R1          | Vertical stress at surface                 | Load center                | MPa    |
| R2          | Horizontal strain at the bottom of layer 1 | ''                         | micron |
| R3          | Vertical strain at the top of layer 2      | ''                           | micron |
| R4          | Vertical strain at the top of layer 3      | ''                           | micron |
|-------------|--------------------------------------------|----------------------------|--------|
| R5          | Vertical stress at surface                 | Edge of load               | MPa    |
| R6          | Horizontal strain at the bottom of layer 1 | ''                           | micron |
| R7          | Vertical strain at the top of layer 2      | ''                           | micron |
| R8          | Vertical strain at the top of layer 3      | ''                           | micron |
<p>
<b>Table 2</b>: Description of response key-points used in basic validation.
</p>

| Software    |   R1    |     R2   |    R3   |    R4   |  R5   |     R6   |  R7     |    R8  |
|-------------|---------|----------|---------|---------|-------|----------|---------|--------|
| BISAR       |   0.7   |  -100.5  |  251.7  |    185  |  0.4  |   -61.9  |  192.2  |  177.5 |
| KENLAYER    |   0.8   |  -100.5  |  251.6  |  185.3  |  0.3  |     -62  |  192.2  |    177 |
| GAMES       |   0.7   |  -100.5  |  251.6  |  185.1  |  0.3  |   -61.9  |  192.2  |  177.5 |
| ALVA (slip) |   0.7   |  -100.4  |  251.6  |  185.1  |  0.3  |   -61.8  |  192.3  |  177.5 |
| ALVA (bond) |   0.7   |  -100.4  |  251.6  |  185.1  |  0.3  |   -61.8  |  192.3  |  177.5 |
<p>
<b>Table 3</b>: Reference pavement system subjeceted to a single wheel load - bonded interfaces.
</p>

| Software    |   R1    |     R2   |    R3   |    R4   |  R5   |     R6   |  R7     |    R8  |
|-------------|---------|----------|---------|---------|-------|----------|---------|--------|
| BISAR       |   0.7   |    -120  |    1    |   217   | 0.4    |   -78   |  -10    |  205 |
| KENLAYER    |   0.7   |    -120  |    1    |   216   |   0    |   -78   |  -10    |  205 |
| GAMES       |   0.7   |    -119  |   11    |   217   |   0    |   -77   |   -1    |  205 |
| ALVA (slip) |   0.7   |  -119.5  |  0.7    | 216.5   | 0.3    | -77.8   | -9.7    |204.7 |

<p>
<b>Table 4</b>: Reference pavement system subjeceted to a single wheel load - unbonded interface between layer 1 and 2, bonded interface between layer 2 and 3.
</p>

| Software    |   R1    |     R2   |    R3   |    R4   |  R5   |     R6   |  R7     |    R8  |
|-------------|---------|----------|---------|---------|-------|----------|---------|--------|
| BISAR       |   0.7   |    N/A   | 186     |  170    |0      |   N/A    |  182    |  177   |
| KENLAYER    |   1.5   |    -85   | 186     |  170    |0      |   -89    |  183    |  177   |
| GAMES       |   0.7   |    -85   | 186     |  170    |0      |   -89    |  183    |  177   |
| ALVA (slip) |   0.7   |  -84.9   | 186     | 169.7   | 0     |  -88.8   | 183.2   | 177.2  |
| ALVA (bond) |   0.7   |  -84.9   | 186     | 169.7   | 0     |  -88.8   | 183.2   | 177.2  |
<p>
<b>Table 5</b>: Reference pavement system subjeceted to a dual wheel load - bonded interfaces.
</p>

| Software    |   R1    |     R2   |    R3   |    R4   |  R5   |     R6   |  R7     |    R8  |
|-------------|---------|----------|---------|---------|-------|----------|---------|--------|
| BISAR       |   0.7   |     N/A  |     9   |   193   | 0     |   N/A    |  -12    |   204 |
| KENLAYER    |   0.7   |    -120  |    -1   |   216   | 0     |   -78    |  -10    |   205 |
| GAMES       |   0.7   |    -101  |    -3   |   194   | 0     |  -106    |    1    |   205 |
| ALVA (slip) |   0.7   |  -102.6  |  -9.2   | 193.4   | 0     |  -107    |-11.5    | 204.1 |

<p>
<b>Table 6</b>: Reference pavement system subjeceted to a dual wheel load - unbonded interface between layer 1 and 2, bonded interface between layer 2 and 3.
</p>

* `ALVA_let_validation1.m` - tests the implementation of the ALVA LET model calculating  vertical stresses and displacements with depth for a half-space subjected to a single circular load. The results are compared to the analytical Boussinesq solution and the computer programme [ELLEA1](https://orbit.dtu.dk/en/publications/ellea1-isotropic-layered-elasticity-in-excel-pavement-analysis-to). 

<i>Note</i>: Minimum two layers is required for analysis of pavement systems in ALVA. For analysis of one-layer / half-space systems: select identical parameters for each layer, as well as <i>"bonded"</i> interface conditions 

<div>
<img src="images/Validation1.png" width="75%">
</div>

* `ALVA_let_validation2.m` - tests the implementation of the ALVA LET model, calculating stresses and displacements with depth for a multilayered pavement subjected to two circular loads utilizing the method proposed. The results are compared to the computer programme [ELLEA1]( https://orbit.dtu.dk/en/publications/ellea1-isotropic-layered-elasticity-in-excel-pavement-analysis-to).
<div>
<img src="images/Validation2.png" width="75%">
</div>

* `ALVA_let_validation3.m` - tests the implementation of the acclerated ALVA LET model, calculating the surface displacements with length for a multilayered pavement subjected to two circular loads utilizing the method proposed in [Andersen et al. (2020)]( https://orbit.dtu.dk/en/publications/efficient-reevaluation-of-surface-displacements-in-a-layered-elas). The results are compared to the computer programme [ELLEA1]( https://orbit.dtu.dk/en/publications/ellea1-isotropic-layered-elasticity-in-excel-pavement-analysis-to).
<div>
<img src="images/Validation3.png" width="75%">
</div>

* `ALVA_let_validation4.m` - tests the implementation of the ALVA LET model calculating shear stresses with depth for a half-space subjected to a single circular load. The results are compared to the computer programme [ELLEA1]( https://orbit.dtu.dk/en/publications/ellea1-isotropic-layered-elasticity-in-excel-pavement-analysis-to).
<div>
<img src="images/Validation4.png" width="75%">
</div>

* `ALVA_let_validation5.m` - tests the implementation of the ALVA LET model calculating shear strain with depth for a multilayered pavement subjected to two circular loads. The results are compared to the computer programme [ELLEA1]( https://orbit.dtu.dk/en/publications/ellea1-isotropic-layered-elasticity-in-excel-pavement-analysis-to).
<div>
<img src="images/Validation6.png" width="75%">
</div>

* `ALVA_visco_validation1.m` - tests the implementation of the ALVA VE model, calculating the displacements for a single evaluation point on the surface of a multilayered pavement considering a single circular load. The results are compared to the computer programme [ELLVA1]( https://orbit.dtu.dk/en/publications/ellva1-isotropic-layered-viscoelasticity-in-excel-moving-load-adv), see details in [Levenberg (2016)]( https://orbit.dtu.dk/en/publications/viscoelastic-pavement-modeling-with-a-spreadsheet).
<div>
<img src="images/Validation5.png" width="75%">
</div>

## User cases
### Case 1 - Backcalculation
* `ALVA_let_backcalculation.m` - tests the implementation of the ALVA LET model for inferring layer moduli, or so called "backcalculation" of layer moduli, based on Falling Weight Deflectometer (FWD) measurements.

#### Input 

| Load     | Radius (mm)    | Pressure, q (MPa) | x-position (mm) |  y-position (mm) |    
|----------|----------------|-------------------|-----------------|------------------|
| Single   |   150       |     0.7           |      0          |   0              |

<p>
<b>Table 7</b>: Falling Weight Deflectometer load. 
</p>

| Geophone                       ---|   G1    |     G2   |    G3    |    G4    |   G5    |     G6   |    G7     |    G8   | G9     |
|-----------------------------------|---------|----------|----------|----------|---------|----------|-----------|---------|--------|
| Location along x-axis (mm):       |   R1    |    200   |    300   |    450   |  600    |    900   |  1200     |   1500  |  1800  |
| Displacement measurement (&mu;m): | 298.9   |    244.2 |    220.07|    175.0 | 138.0   |   97.3   |   72.2    |     57.4|    47.6|

<p>
<b>Table 8</b>: Sensor location and displacement measurements.
</p>

## Installation
* Download the package on your PC. 
* Open MATLAB
* Go to the directory 'ALVA'
* add the different directories of the ALVA on your MATLAB path — Now you are ready to run the validation examples provided and generate your own analysis. 
* ALAVA is compatible with [OCTAVE](https://www.gnu.org/software/octave/index)
