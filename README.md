Skar, Asmus, Andersen, Sebastian and Julius Nielsen (2020): Adaptive Layered Viscoelastic Analysis (ALVA). Technical University of Denmark. Software. https://doi.org/10.11583/DTU.12387305

# ALVA
Adaptive Layered Viscoelastic Analysis (ALVA) is a MATLAB-package for pavement modeling. The aim with the software is to equip the civil engineering community with an advanced pavement modeling tool and computer package that is highly adaptive, transparent and open-access, capable of supporting current and future pavement evaluation needs. 

<div>
<img src="images/Roadmap_ALVA.png" width="90%">
</div>

The components of the ALVA software (shown in the Figure) is briefly described below:
* `main.m` - main script for defining vehicle loading conditions, pavement structure geometry and material properties, numerical parameters (for balancing accuracy and efficiency), analysis type and evaluation points (i.e., location of the output response), and post-processing of the results.
* `init_LET.m` - initialize the analysis of a layered elastic half space selecting between: (a) `let_response_full.m` - evaluates the response (displacements, stresses and strains) of a layered elastic half-space model or (b) `let_response_polfit.m` - evaluates the displacements of a layered elastic half space model at the surface only. 
* `init_LET.m` also works as an interface application for simulating moving loads in the *Viscoelastic module*; based on the asphalt material properties from `VE_moduli.m` the module is engaged to calculate the elastic responses in `VE_response.m`, which is subsequently used to predict the linear viscoelastic response of a moving vehicle in `VE_simulation.m`.
* The main support functions in `init_LET.m` are (a) `arb_funct.m` - evaluates the coefficients of integration that embody the layered system properties and connectivity and (b) `numint_coeff.m` – organizes the integration points and weights for numerical integration for the user defined number of Bessel roots from `besselroots.m` and gauss points and weights between Bessel roots from `lookup_gauss.m`.

* Additional support functions for efficient surface displacement calculations in `let_response_polfit.m` is `arb_func_polfit.m` – interpolates intermediate coefficients of integration, `polfit_int.m` – evaluates coefficient proportional integrals, and `polfit_abc.m` – evaluates coefficients a, b and c for the polynomial.
Examples on ‘main.m’ scripts can be found in the _../ALVA/examples_ folder and all underlying components (i.e., source files) described above in the _../ALVA/basic_ folder.

## Main input parameters
The core algorithm behind this package is based on LET, i.e., the classic formulation for an N-layered half-space, shown below.

<div>
<img src="images/N_layer.png" width="90%">
</div>

In this model all layers are assumed linear elastic, isotropic, homogeneous, fully bonded, and weightless. The model inputs include Young’s modulus E<sub>n</sub>, Poisson’s ratio &Nu;<sub>n</sub>, and layer thickness t<sub>n</sub> (where n denotes the layer number). This model is engaged to calculate the response at any point, A<sub>j</sub>, of interest and for a given set of uniformly distributed circular loadings with load radius, a, and pressure q). An overview of LET model assumptions and solution procedure is given in [Khazanovich and Wang (2007)](https://journals.sagepub.com/doi/abs/10.3141/2037-06)

The viscoelastic response is approximated based on the LET calculations utilizing the methodology and load scheme suggested by [Levenberg (2016)](https://orbit.dtu.dk/en/publications/viscoelastic-pavement-modeling-with-a-spreadsheet). Viscoelastic layers are associated with a creep compliance and model parameters D<sub>0</sub> and D<sub>&infin;</sub>, the short and long time compliances (respectively), and shape parameters &tau;<sub>D</sub> and n<sub>D</sub>, controlling the transition between D<sub>0</sub> and D<sub>&infin;</sub>. 

<div>
<img src="images/VE_mesh.png" width="90%">
</div>

The load moves in a straight line from x=-x<sub>0</sub> (Start) to x=x<sub>0</sub> (End). The travel path is decomposed into N intervals (i=1,…,N), each &Delta;x long. The point of response evaluation A<sub>j</sub> is indicated in the Figure; this point is located near the middle of the travel path (i.e., x-coordinate of zero), at y-coordinate y<sub>0</sub> and depth z<sub>0</sub> below the surface. 

## Validation examples
ALVA comes with seven validation examples (i.e., `main.m` scripts) comparing ALVA to existing codes and analytical formulations. 
The example files can be found in the in the _../ALVA/examples_ folder and the results obtained with independent codes can be found in the found in the _../ALVA/validation_ folder.  
The following examples are included in this package:

* 'ALVA_bonding_validation1.m' - tests the implementation of the interface spring model, as well as compares the ALVA model with a range of commonly used softaware at critical positions within the pavament system published by the [European Commission](https://trimis.ec.europa.eu/project/advanced-models-analytical-design-european-pavement-structures).


| Layer    | Thickness (mm) | Youngs moduli (MPa) | Poisson's ratio  |    
|----------|----------------|---------------------|------------------|
| 1        |   260          |     5000            |   0.35           |
| 2        |   500          |      200            |   0.40           |
| 3        |   &infin;      |       50            |   0.45           |

<p>
<b>Table 1</b>: Reference pavement system
</p>

| Load     | Radius (mm)    | Pressure, q (MPa) | x-position (mm) |  y-position (mm) |    
|----------|----------------|-------------------|-----------------|------------------|
| Single   |   150.8        |     0.7           |      0          |   0              |
| Dual (1) |   106.6        |     0.7           |   -170          |   0              |
| Dual (2) |   106.6        |     0.7           |    170          |   0              |


<p>
<b>Table 2</b>: Load cases 
</p>


<p>
| Response, R | Description                                        |
|-------------|----------------------------------------------------|
| R1          | Vertical stress surface at center of load          | 
| R2          | Horizontal strain bottom layer 1 at center of load |
| R3          | Vertical strain top layer 2 at center of load      | 
| R4          | Vertical strain top layer 3 at center of load      | 
| R5          | Vertical stress surface at edge of load            | 
| R6          | Horizontal strain bottom layer 1 at edge of load   | 
| R7          | Vertical strain top layer 2 at edge of load        | 
| R8          | Vertical strain top layer 3 at edge of load        |
</p>
<p>
<b>Table 2</b>: Description of response key-points for validations
</p>

| Software    |   R1    |     R2   |    R3   |    R4   |  R5   |     R6   |  R7     |    R8  |
|-------------|---------|----------|---------|---------|-------|----------|---------|--------|
| BISAR       |   0.7   |  -100.5  |  251.7  |    185  |  0.4  |   -61.9  |  192.2  |  177.5 |
| KENLAYER    |   0.8   |  -100.5  |  251.6  |  185.3  |  0.3  |     -62  |  192.2  |    177 |
| GAMES       |   0.7   |  -100.5  |  251.6  |  185.1  |  0.3  |   -61.9  |  192.2  |  177.5 |
| ALVA (slip) |   0.7   |  -100.4  |  251.6  |  185.1  |  0.3  |   -61.8  |  192.3  |  177.5 |
| ALVA (bond) |   0.7   |  -100.4  |  251.6  |  185.1  |  0.3  |   -61.8  |  192.3  |  177.5 |
<p>
<b>Table 3</b>: Reference pavement system subjeceted to a single wheel load - bonded interfaces
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

## Installation
* Download the package on your PC. 
* Open MATLAB
* Go to the directory 'ALVA'
* add the different directories of the ALVA on your MATLAB path — Now you are ready to run the validation examples provided and generate your own analysis. 
