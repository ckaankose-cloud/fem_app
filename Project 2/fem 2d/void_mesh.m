  function [NL, EL] = void_mesh(d1, d2, p, m, R, element_type)

q = [0 0; d1 0; 0 d2; d1 d2];

PD = 2;
NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1);
NoE = 4*p*m;
NPE = 4;


%%% Nodes %%%

NL = zeros(NoN, PD);

a = d1/p;
b = d2/p;

%%%% Region 1  %%%

coor11 = zeros((p+1)*(m+1), PD);

for i = 1:p+1
    
    coor11(i,1) = q(1,1) + (i-1)*a;
    coor11(i,2) = q(1,2);
    
end

for i = 1:p+1
    coor11(m*(p+1)+i,1) = R*cos((5*pi/4) + (i-1)*(pi/2)/p) + (d1/2);
    coor11(m*(p+1)+i,2) = R*sin((5*pi/4) + (i-1)*(pi/2)/p) + (d2/2);
end


for i = 1:m-1
    
    for j = 1:p+1
        
        dx = (coor11(m*(p+1)+j,1) - coor11(j,1))/m;
        dy = (coor11(m*(p+1)+j,2) - coor11(j,2))/m;
        
        coor11(i*(p+1)+j,1) = coor11(((i-1)*(p+1))+j,1) + dx;
        coor11(i*(p+1)+j,2) = coor11(((i-1)*(p+1))+j,2) + dy;
        
    end
    
end


%%%% Region 2  %%%


coor22 = zeros((p+1)*(m+1), PD);


for i = 1:p+1
    
    coor22(i,1) = q(3,1) + (i-1)*a;
    coor22(i,2) = q(3,2);
    
end

for i = 1:p+1
    coor22(m*(p+1)+i,1) = R*cos((3*pi/4) - (i-1)*(pi/2)/p) + (d1/2);
    coor22(m*(p+1)+i,2) = R*sin((3*pi/4) - (i-1)*(pi/2)/p) + (d2/2);
end

for i = 1:m-1
    
    for j = 1:p+1
        
        dx = (coor22(m*(p+1)+j,1) - coor22(j,1))/m;
        dy = (coor22(m*(p+1)+j,2) - coor22(j,2))/m;
        
        coor22(i*(p+1)+j,1) = coor22(((i-1)*(p+1))+j,1) + dx;
        coor22(i*(p+1)+j,2) = coor22(((i-1)*(p+1))+j,2) + dy;
        
    end
    
end
   

%%%  Region 3  %%%%

coor33 = zeros((p-1)*(m+1), PD);

for i = 1:p-1
    
    coor33(i,1) = q(1,1);
    coor33(i,2) = q(1,2) + (i)*b;
    
end


for i = 1:p-1
    coor33(m*(p-1)+i,1) = R*cos((5*pi/4) - (i)*(pi/2)/p) + (d1/2);
    coor33(m*(p-1)+i,2) = R*sin((5*pi/4) - (i)*(pi/2)/p) + (d2/2);
end

for i = 1:m-1
    
    for j = 1:p-1
        
        dx = (coor33(m*(p-1)+j,1) - coor33(j,1))/m;
        dy = (coor33(m*(p-1)+j,2) - coor33(j,2))/m;
        
        coor33(i*(p-1)+j,1) = coor33(((i-1)*(p-1))+j,1) + dx;
        coor33(i*(p-1)+j,2) = coor33(((i-1)*(p-1))+j,2) + dy;
        
    end
end


%%%%  Region 4  %%%%

coor44 = zeros((p-1)*(m+1), PD);

for i = 1:p-1
    
    coor44(i,1) = q(2,1);
    coor44(i,2) = q(2,2) + (i)*b;
    
end

for i = 1:p-1
    
    coor44(m*(p-1)+i,1) = R*cos((7*pi/4) + (i)*(pi/2)/p) + (d1/2);
    coor44(m*(p-1)+i,2) = R*sin((7*pi/4) + (i)*(pi/2)/p) + (d2/2);
    
end

for i = 1:m-1
    
    for j = 1:p-1
        
        dx = (coor44(m*(p-1)+j,1) - coor44(j,1))/m;
        dy = (coor44(m*(p-1)+j,2) - coor44(j,2))/m;
        
        coor44(i*(p-1)+j,1) = coor44(((i-1)*(p-1))+j,1) + dx;
        coor44(i*(p-1)+j,2) = coor44(((i-1)*(p-1))+j,2) + dy;
        
    end
    
end


for i = 1:m+1 
    
    NL((i-1)*4*p+1:(i)*4*p,:) = [coor11((i-1)*(p+1)+1:(i)*(p+1),:);
                                 coor44((i-1)*(p-1)+1:(i)*(p-1),:);
                                 flipud(coor22((i-1)*(p+1)+1:(i)*(p+1),:));
                                 flipud(coor33((i-1)*(p-1)+1:(i)*(p-1),:))];
    
end

%buraya kadar nodeları sıraladı

%%%%   Elements   %%%%


EL = zeros(NoE, NPE);

for i = 1:m 
    
    for j = 1:4*p
        
        if j == 1
            
            EL((i-1)*(4*p)+j,1) = (i-1)*(4*p) + 1;
            EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) + 1;
            EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p;
            EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1;
            
        elseif j == 4*p
            
            EL((i-1)*(4*p)+j,1) = (i)*(4*p);
            EL((i-1)*(4*p)+j,2) = (i-1)*(4*p) + 1;
            EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,1) + 1;
            EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p;
            
        else
            
            EL((i-1)*(4*p)+j,1) = EL((i-1)*(4*p)+j-1,2);
            EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j-1,3);
            EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1;
            EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) + 1;
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isequal(element_type,'D2TR3N')
    
    NPE_new = 3;
    EL_new = zeros(2*NoE, NPE_new);
    
    for i = 1:NoE
        
        EL_new(2*(i-1)+1,1) = EL(i,1); 
        EL_new(2*(i-1)+1,2) = EL(i,2);
        EL_new(2*(i-1)+1,3) = EL(i,3);
        
        EL_new(2*(i-1)+2,1) = EL(i,1); 
        EL_new(2*(i-1)+2,2) = EL(i,3);
        EL_new(2*(i-1)+2,3) = EL(i,4);
        
    end
    
    EL = EL_new;
    
end


end