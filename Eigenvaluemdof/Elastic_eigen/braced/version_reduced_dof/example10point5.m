function [nodes,StoreyNodes,elems,bcs,bcs_unbraced,loads]=example10point5
    

nodes=[0.0 	0.0;
        0   288;
        240	288;
        720 288;
        720  0;
        240  0];
    
StoreyNodes=[1 6 5;
             2 3 4];    

                          %(Node1 Node2), E, I, A,Fy,Zx
elems= [1 2 29000 248 13.3 36 54.9;
        6 3 29000 248 13.3 36 54.9;
        5 4 29000 248 13.3 36 54.9;
        2 3 29000 2850 24.8 36 244;
        3 4 29000 2850 24.8 36 244];

                       % (Node  dof  disp)
bcs=[1 1 0;
     1 2 0;
     1 3 0;
     6 1 0;
     6 2 0;
     6 3 0;
     5 1 0;
     5 2 0;
     5 3 0];
% (Node Storey dir  val)
 bcs_unbraced=[1 1 1 0;
               1 1 2 0;
               1 1 3 0;
               6 1 1 0;
               6 1 2 0;
               6 1 3 0;
               5 1 1 0;
               5 1 2 0;
               5 1 3 0];

                  % (Node  dof  load)
loads=[2 2 -1;
       3 2 -1;
       4 2 -1];

end