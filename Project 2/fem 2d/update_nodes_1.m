function [ENL] = update_nodes_1(ENL,Uu,NL,Fu)

NoN = size(NL,1);
PD = size(NL,2);
DOFs = 0; 
DOCs = 0;
    
for i=1:NoN
    for j=1:PD
        if ENL(i,PD+j)==1 %if bc is neumann,force is described
            DOFs=DOFs+1;
            ENL(i,4*PD+j)=Uu(DOFs);% fill in the unknown displacement 
        end
    end
end

for i=1:NoN
    for j= 1:PD
        if ENL(i,PD+j)==-1 % if bc is dirichlet, displacement described 
            
            DOCs=DOCs+1;
            ENL(i,5*PD+j)=Fu(DOCs);%fill in the unknown forces 
        end
    end
end

end