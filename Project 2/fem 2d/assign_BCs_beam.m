function [ENL , DOFs, DOCs] = assign_BCs_beam(NL, BC_type, magnitude, d1, d2)

NoN = size(NL,1); % Number of Nodes
PD = size(NL,2); %Problem Dimension

ENL = zeros(NoN,PD*6);

ENL(:,1:PD) = NL; %Coordinates

switch BC_type
    case 'y_loading'

        for i = 1 : NoN
            
            if ENL(i,1) == 0
                
                ENL(i,PD+1) = -1;
                ENL(i,PD+2) = -1;
                
                ENL(i,4*PD+1) = 0;
                ENL(i,4*PD+2) = 0;
                
            elseif ENL(i,1) == d1
                
                ENL(i,PD+1) = -1;
                ENL(i,PD+2) = -1;
                
                ENL(i,4*PD+1) = 0;
                ENL(i,4*PD+2) = -magnitude;
            
            else
                ENL(i,PD+1) = 1;
                ENL(i,PD+2) = 1;
                
                ENL(i,4*PD+1) = 0;
                ENL(i,4*PD+2) = 0;
            end
        end
    
        
   
    case 'x_loading'
        
        for i = 1 : NoN
            
            if ENL(i,1) == 0
                
                ENL(i,PD+1) = -1;
                ENL(i,PD+2) = -1;
                
                ENL(i,4*PD+1) = 0;
                ENL(i,4*PD+2) = 0;
                
            elseif ENL(i,2) == d2
                
                ENL(i,PD+1) = -1;
                ENL(i,PD+2) = -1;
                
                ENL(i,4*PD+1) = magnitude;
                ENL(i,4*PD+2) = 0;
                
            else
                ENL(i,PD+1) = 1;
                ENL(i,PD+2) = 1;
                
                ENL(i,4*PD+1) = 0;
                ENL(i,4*PD+2) = 0;
            end
        end

end

DOFs = 0;
DOCs = 0;
       
for i = 1:NoN
    
    for j = 1:PD
        
        if (ENL(i,(PD+j)) == -1)
            
            DOCs = DOCs - 1;
            ENL(i,2*PD+j) = DOCs;
            
        else
            
            DOFs = DOFs + 1;
            ENL(i,2*PD+j) = DOFs;
        end
    end
end

for i = 1:NoN
    
    for j = 1:PD
    
        if (ENL(i,(2*PD+j))<0)
            
            ENL(i,PD*3+j) = DOFs + abs(ENL(i,(2*PD+j)));
            
        else
            
            ENL(i,PD*3+j) = ENL(i,(2*PD+j));
            
        end
    end
end

DOCs = abs(DOCs);

end
    
      
   