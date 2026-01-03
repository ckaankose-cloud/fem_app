function [K] = assemble_stiffness(ENL, EL, NL,E,nu)

NoE = size(EL,1);
NPE = size(EL,2); 

NoN = size(NL,1);
PD = size(NL,2);

K = zeros(NoN*PD, NoN*PD);

if NPE == 3
    GPE = 1;
elseif NPE == 4
    GPE = 4; 
elseif NPE == 6
    GPE = 3;
elseif NPE == 8
    GPE = 4;
elseif NPE ==9
    GPE = 9;  
end

for i = 1:NoE
    
    nl = EL(i,1:NPE); 
    
    x = zeros(NPE, PD);
    
    for j = 1:NPE
        
        x(j,:) = ENL(nl(j),1:PD);

    end

k = Stiffness(x, GPE, E, nu);  % Stiffness calculates element stiffness
    
    for r = 1:NPE %d1
        
        for p = 1:PD
           
            for q = 1:NPE %d2
                
                for s = 1:PD
                    row = ENL(nl(r), p+3*PD); %Orders wrt the global DOFs
                    column = ENL(nl(q),s+3*PD);
                    
                    value = k((r-1)*PD+p, (q-1)*PD+s);
                    K(row,column) = K(row,column) + value;
                end
            end
        end
    end
end
                