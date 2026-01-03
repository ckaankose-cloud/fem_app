function [Up] = assemble_displacement_1(ENL,NL)
    
NoN = size(NL,1); %12
PD  = size(NL,2); %2 

DOC=0;

for i=1:NoN
    for j =1:PD
        if ENL(i, PD+j)==-1  
            DOC=DOC+1;
            Up(DOC,1) = ENL(i,4*PD+j);
        end
    end
end

end