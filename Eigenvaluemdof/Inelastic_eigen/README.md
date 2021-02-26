This Module is written to perform INELASTIC Eigenvalue Analysis using Cubic Hermitian Elements and Power Method. It contains the main script Inelastic_eigen.m that calls upon
the lamdainelastic.m function to obtain a linearized estimate of inelastic buckling load. Regular Falsi iterations are then performed to obtain the final converged buckling load. 
The script allows detection of cases where elastic buckling governs i.e. axial stresses are below proportional stress (Fy-Fresidual). 
