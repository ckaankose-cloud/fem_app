function [K] = assemble_stiffness_1(ENL, EL , NL, E, A) %takes as input 

NoE=size(EL,1); %number of rows , number of elements
NPE=size(EL,2); %number of columns, nodes per element

NoN=size(NL,1);%number of nodes
PD=size(NL,2);%problem dimensions
    
K=zeros(NoN*PD,NoN*PD);%NoNxPD

for i=1:NoE
    nl=EL(i,1:NPE);%node list of each element for i=1 1 2, for i=2 2 3,for i=3 1 3.
    k = element_stiffness_1(nl,ENL,E,A);%takes input nodelist of each element 
    for r=1:NPE %NPE 
        for p=1:PD %PD 
            %%%%%%%%%%%%%%
            for q=1:NPE %d2
                for s=1:PD %PD
                    row= ENL(nl(r),p+3*PD);%nl(1),7 =4
                    column= ENL(nl(q),s+3*PD);%nl(1)
                    
                    value= k((r-1)*PD+p, (q-1)*PD+s);
                    K(row,column)= K(row,column) + value;
                end
            end
        end
    end
end