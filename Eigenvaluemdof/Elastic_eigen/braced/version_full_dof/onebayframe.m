function [nodes,elems,bcs,loads]=onebayframe
    

nodes=[0.0 	0.0;
        0   96;
        48	96];
        
        

%     (N1N2) E,   I,    A,   Fy, Zx
elems= [1 2 29000 248  13.3 36 54.9;
        2 3 29000 248  13.3 36 54.9];
        
        

                       % (Node  dof  disp)
bcs=[1 1 0;
     1 2 0;          
     3 1 0;
     3 3 0];
          
     
     
%{
                  % (Node  dof  load)
loads=[2 1 0.1;
       3 2 -1;
       4 2 -0.5];
%}
  loads=[2 2 -1];
         
            
     
     
end