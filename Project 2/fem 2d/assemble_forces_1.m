function [Fp]= assemble_forces_1(ENL,NL) %ENL and L input

NoN = size(NL,1);
PD  = size(NL,2);

DOF = 0;

for i=1:NoN   
    for j = 1:PD
        if ENL(i,PD+j) == 1 % if bc is neumann

            DOF= DOF+1; %inrcease dof
            Fp(DOF,1) = ENL(i,5*PD+j); %add the forces to Fp
            
        end
    end
end

end