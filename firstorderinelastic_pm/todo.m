TO DO: 

1a:Include Natural Deformation to euler/rk2 (elastic) and compare against NRM
1b: Include Natural Deformation state determination to euler/rk2 inelastic code and compare 
    against standard force recovery method. 

2:FOR NRM, consider the following iterative state determination schemes for 
force recovery:
version1: use nodexy and Pj from last iterate state to form kt-local and 
          use standard deformation vector to compute dfint. 
version2: use nodexy and Pj from last iterate state to form kt-local and 
          use natural deformation vector to compute dfint.
version3: use midpoint Pj and nodexy (computed for the given delta increment/2)
          to form secant stiffness matrix and use that with/without natural
          deformations vector to compute dfint.
version4: use secant stiffness matrix (co-rotational formulation given in Chen
          and Sohal Book) to evaluate corotational P-Ma-Mb state through iterations
          and then use that for Finternal determination. 

Compare MNRM and NRM performance in terms of iterations, step length, convergence 
characteristics etc.

3:Include provision for AUTOMATED STEP LENGTH and improved convergence criteria in NRM.

Compare Cubic Hermitian and Stabililty Function Kg matrics. 
Write up codes for Elastic/Inelastic Eigenvalue Critical Load Analysis

4a:Try to incorporate mixed NRM+rk2 approach for running second order inelastic analysis 
of frame structures. (one approach would be to run standard NRM until first hinge formation
and then switch to typical RK-2 scheme for rest of the analysis.) Another approach would be
to run full NRM all the way. Lets see! 

4b:Include effect of member/node discretization and imperfection in the current codes. 

4c: Compare effect of  proportional and non-proportional loading on nonlinear analysis results.

5:Upgrade to 3D for capturing mixed bi-axial flexural/torsional deformation modes. 

PAPER WRITEUP!!!
THEEND





Comparison of Alignment Charts Method to Eigenvalue Analysis??
Compare Direct Design Philosophy AISC vs Approximate Second Order Amplification Methods.

Look into Distributed Plastic Zone method and residual stress inclusion in non-linear analysis. 

5z:consider using exact P-M interaction curve and compare with Ziemian convex yield surface approximation 10.18  (Orbinson McGuire Ref 10.10)
5b:variable load step during analysis. 
5c1:inclusion of tangent modulus to account for residual stresses and partial yielding(flexure) (2nd order inelastic analysis): for axial and flexural deformations
5c2:look into Et,Etm option in MASTAN and quasi hinge modelling of Ziemian chapter 10.
distributed plasticitu Alemdar Bulent .
5d:inclusion of element interior nodes to account for p-small delta
5e:inclusion of member imperfections 
5f:inclusion of shear deformations
5g:inclusion of partially rigid connections in non-linear analysis.
5h:inclusion of panel zone deformations
5i:Effect of actual end-offsets of member lengths
5j:Effect of proportional and non-proportional loading
5k:Plotting post-peak response curves for 2ndorder-inelastic response(only applicable if K-structure at peak load is not positive definite and NOT ill-conditioned. 
5l:Accounting for settlement induced nonlinear response

6a:fibre hinge mixed finite element approach to model frame members

7a:Cyclic/Dynamic Inelasticity: kinematic/isotropic hardening
   Newmark algorithms and energy conserving largrange multipliers
   
8a:OpenSees Modelling of Steel Frames. Foray into Braced Frames Response.   


