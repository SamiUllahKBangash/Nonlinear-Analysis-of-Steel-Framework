function [nodes,StoreyNodes,elems,bcs,bcs_unbraced,loads]=Fullonebayframe
    

nodes=[0.0 	0.0;
        0   96;
        0	192;
        96  192;
        96  96;
        96   0];
        
        
StoreyNodes=[1 6;
             2 5;
             3 4];
    
%     (N1N2) E,   I,    A,   Fy, Zx
elems= [1 2 29000 248  13.3  36 54.9;
        2 3 29000 248  13.3  36 54.9;
        3 4 29000 248  13.3  36 54.9;
        4 5 29000 248  13.3  36 54.9;
        5 6 29000 248  13.3  36 54.9;
        2 5 29000 248  13.3  36 54.9];
        
        

% (Node dir  val)
bcs=[1 1 0;
     1 2 0;
     1 3 0;
     6 1 0;
     6 2 0
     6 3 0];
          
% (Node Storey dir  val)
 bcs_unbraced=[1 1 1 0;
               1 1 2 0;
               1 1 3 0;
               6 1 1 0;
               6 1 2 0;
               6 1 3 0];              
                  
     
%{
                  % (Node  dof  load)
loads=[2 1 0.1;
       3 2 -1;
       4 2 -0.5];
%}
  loads=[3 2 -1;
         4 2 -1];
      
         
            
     
     
end