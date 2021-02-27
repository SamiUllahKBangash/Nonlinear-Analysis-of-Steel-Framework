This module performs first order inelastic analysis using event to event plastic hinge procedure. Bilinear elastic-plastic stress strain response is assumed
with user-defined hardening parameter. Axial-Flexural plastic interaction is not accounted for i.e. (plasticity is only assumed on internal moments at element ends') Co-rotational
and Coordinate Transformations are used to obtain Global-Local Stiffness Matrices based on the work by Chen et al. Element State Determination is used to obtain successive Load
Scale Factors to jump to the next plastic hinge event(s). 

Instructions:
Upon running Final_1storderinelastic_cleaned.m, a pop up window will appear from where one can choose ChenBookExample.txt input file for extracting frame geometry. You may need to reduce plot variable range by 1 to see plastic hinge formation more clearly. 

Acknowledgements: Asghar Jadoon, NICE, NUST
