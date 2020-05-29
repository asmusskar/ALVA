---
title: 'ALVA: An adaptive MATLAB package for layered viscoelastic elastic analysis'
tags:
  - MATLAB
  - Layered Elastic Theory
  - Linear Viscoelasticity
  - Pavement analysis
authors:
  - name: Asmus Skar
    orcid: 0000-0003-3176-791X
    affiliation: 1
  - name: Sebastian Andersen
    orcid: 0000-0002-6397-7914
    affiliation: 1
affiliations:
 - name: Department of Civil Engineering, Technical University of Denmark, 2800 Kgs. Lyngby, Denmark
   index: 1
date: 29 May 2020
bibliography: paper.bib
---

# Summary
A key component in design and analysis of asphalt pavements is the response model, used to calculate the stresses, strains and displacements in the structure, when subjected to mechanical loading. For many years, response models formulated within the theoretical continuum-mechanics framework of a stratified half-space have been successfully used [@khazanovich:2007a]. In this framework a pavement system is viewed as a collection of layers, each of finite thickness, resting on a semi-infinite medium, also referred to as the Layered Elastic Theory (LET) [@Burmister:1945a].  Several computer programs for the analysis of pavement systems were developed on the basis of LET. However, new design guidelines and in-situ evaluation of mechanical pavement properties entail re-execution of the underlying model many times over, making the response analysis extremely time consuming. This has resulted in development of numerical acceleration techniques for improved computational efficiency [@Erlingsson:2013a;@khazanovich:2007a] implemented in the computer package ELLEA [@Levenberg:2016a]. Moreover, the engineering community is placing increased emphasis on linear viscoelastic characterization of asphalt materials. It is therefore anticipated, that routine pavement-related calculations would soon evolve to include time dependent layer properties and moving loads as supported by the computer packages ViscoRoute [@Chabot:2010a] and ELLVA [@Levenberg:2016b]. Other relevant features for pavement evaluation is the capability of describing both perfect and deficient bonding along the interfaces of neighboring layers [@Maina:2007a].

Although, pavement analysis tools have been developed to address the limitations of conventional LET analysis, there are currently no computer package that supports all aforementioned needs in a unified manner. Moreover, existing software are often limited to standard design problems (e.g., only applicable to a limited number of layers and loads) and the control over the program features is restricted (e.g., source codes not available, implementation not documented, numerical parameters fixed). Thus, the aim with the software presented herein is to equip the civil engineering community with an advanced pavement modeling tool and computer package that is highly adaptive, transparent and open-access, capable of supporting current and future pavement evaluation needs. 

To achieve this, we have implemented some of the most promising numerical techniques currently available for effectively solving elastic and linear viscoelastic layered problems in a MATLAB computer package. The core algorithm behind this package is based on LET, including numerical features for improved code performance, e.g., acceleration techniques to speed computational time for analysis of surface displacements [@Andersen:2020a] and full response analysis [@khazanovich:2007a; @Levenberg:2016a] and extrapolation technique to improve convergence [@Erlingsson:2013a]. The viscoelastic response is approximated based on the LET calculations utilizing the methodology proposed by [@Levenberg:2016c]. In addition ALVA allows users to define horizontally oriented springs operating at the top or bottom layer interfaces to model imperfect interface conditions [@Levenberg:2020a].

The final computer package Adaptive Layered Viscoelastic Analysis (ALVA) offers a near real-time solution for the history of stress, strain, and displacement inside the system at any point of interest resulting from a moving load. The generic computational scheme proposed enables users to control any feature of the program, including geometrical and material properties, layer interface bonding, loading conditions and numerical parameters. Thus, ALVA supports advanced analysis, e.g., detailed analysis of surface tire interaction and unconventional axle loads, and enables users to balance numerical efficiency and accuracy. ALVA can be used as is [@Andersen:2020a;@Skar:2020a] or serve as computational kernel in new software; supporting future advances in pavement analysis software and development of model-guided data interpretation schemes.
Future developments in this package may include mathematical formulations to describe more complex phenomena, e.g., *(i)* fragmented layer conditions (multi-cracked layers); *(ii)* both isotropic or transversely isotropic properties; *(iii)* both vertical and horizontal surface loads. 

# References
