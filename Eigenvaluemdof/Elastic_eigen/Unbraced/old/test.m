[nodes,StoreyNodes,elems,bcs,bcs_unbraced,loads]=Fullonebayframe;
Nel=size(elems,1);          %No of frame elements
Nnodes=size(nodes,1);       %No of Nodes in the structure

alldofs=1:(Nnodes + size(StoreyNodes,1));
dofspec=[]; %constrained DOFs Array
for ii=1:size(bcs_unbraced,1)
    if bcs_unbraced(ii,3)==3
        thisdof=bcs(ii,1);
        dofspec=[dofspec thisdof];
    elseif bcs_unbraced(ii,3)==1
        thisdof=Nnodes + bcs(ii,2);
        if  isempty(dofspec==thisdof)
            dofspec=[dofspec thisdof];
        end
    end
         %u(thisdof)=bcs(ii,3);    %currently,the code doesn't account for nonlinear node settlement cases.
    %du(thisdof)=bcs(ii,3);
end
doffree=alldofs;
doffree(dofspec)=[]; %alldofs= (doffree + dofspec)

