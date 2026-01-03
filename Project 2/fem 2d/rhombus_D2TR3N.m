function [NL, EL] = rhombus_D2TR3N(d1, d2, p, m, e, inclusion)

q = [0 0; d1 0; 0 d2; d1 d2];%corners of the domain

PD = 2;

NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1); %regions 1-2 + regions 3-4

NoE = 4*p*m;

NPE= 4; %quadrilateral

%%% nodes %%%

NL = zeros(NoN, PD);%number of nodes initialized

a = d1/p;%increment in x
b = d2/p;%increment in y
%region 1and 2 include diagonal ,3  and 4 exclude diagonal

%%%% Region 1 bottom %%% from left to right

coor11 = zeros((p+1)*(m+1), PD);    %small node list of region1 
    
for i = 1:p+1 %nodes along column along d1

    coor11(i,1) = q(1,1) + (i-1)*a; %x coordinate of first corner +increment
    coor11(i,2) = q(1,2); %y coordinate of first corner +increment
end

for i = 1:p+1 % nodes along column  along void curve
    coor11(m*(p+1)+i,1) = (d1/2-e/(2*sqrt(2))) + (i-1)*e/(p*sqrt(2)) ;%x component + increment
    if mod(p,2)==0
        mult = (p+1)/2;
    else
        mult = (p+1)/2 + 1;
    end

    if i<mult
        coor11(m*(p+1)+i,2) = (d2/2-e/(2*sqrt(2))) - (i-1)*e/(p*sqrt(2));
    elseif i == (p+2)/2
        coor11(m*(p+1)+i,2) = d2/2 - e/sqrt(2);
    else
        coor11(m*(p+1)+i,2) = (d2/2-e/(2*sqrt(2))) + (i-1-(p))*e/(p*sqrt(2));
    end
%         coor11(m*(p+1)+i,2) = (d2/2-e/(2*sqrt(2))) ;%
end

%generation of the nodes in remaining region%
for i = 1:m-1 %along rows
    
    for j = 1:p+1 %along columns
        
        dx = (coor11(m*(p+1)+j,1) - coor11(j,1))/m;%increment x over m increments
        dy = (coor11(m*(p+1)+j,2) - coor11(j,2))/m;%increment y
        
        coor11(i*(p+1)+j,1) = coor11(((i-1)*(p+1))+j,1) + dx;
        coor11(i*(p+1)+j,2) = coor11(((i-1)*(p+1))+j,2) + dy;
    end
end


%%%% Region 2 top  %%% from left to right


coor22 = zeros((p+1)*(m+1), PD);

for i = 1:p+1
    
    coor22(i,1) = q(3,1) + (i-1)*a; %x coord + increment
    coor22(i,2) = q(3,2); 

end

for i = 1:p+1

    coor22(m*(p+1)+i,1) = ((d1+d2)/4+e/(2*sqrt(2))) + (i-1-p)*e/(p*sqrt(2)) ;%x component + increment
    if mod(p,2)==0
        mult = (p+1)/2;
    else
        mult = (p+1)/2 + 1;
    end
    if i<mult

        coor22(m*(p+1)+i,2) = (d2)/2 + (e/(2*sqrt(2))) + (i-1)*e/(p*sqrt(2));
    elseif i == (p+2)/2
        coor22(m*(p+1)+i,2) = d2/2 + e/sqrt(2);
    else
        coor22(m*(p+1)+i,2) = (d2)/2 + (e/(2*sqrt(2))) - (i-1-(p))*e/(p*sqrt(2));
    end

    %coor22(m*(p+1)+i,1) = ((d1+d2)/4+e/(2*sqrt(2)));
    %coor22(m*(p+1)+i,2) = (d2)/2 + (e/(2*sqrt(2)));
end

for i = 1:m-1
    
    for j = 1:p+1
        
        dx = (coor22(m*(p+1)+j,1) - coor22(j,1))/m;
        dy = (coor22(m*(p+1)+j,2) - coor22(j,2))/m;
        
        coor22(i*(p+1)+j,1) = coor22(((i-1)*(p+1))+j,1) + dx;
        coor22(i*(p+1)+j,2) = coor22(((i-1)*(p+1))+j,2) + dy;
        
    end
    
end
   

%%%  Region 3 left %%%% from bottom to top 

coor33 = zeros((p-1)*(m+1), PD);

for i = 1:p-1
    
    coor33(i,1) = q(1,1);
    coor33(i,2) = q(1,2) + (i)*b;%skip on the node diagonal
    
end

for i = 1:p-1 %along the curve
    
    if mod(p,2)==0
        mult = (p+1)/2;
    else
        mult = (p+1)/2 + 1;
    end

    if i<mult - 1
        coor33(m*(p-1)+i,1) = (d1+d2)/4 - e/(2*sqrt(2)) - (i)*e/(p*sqrt(2));
    elseif i == (p+2)/2 - 1
        coor33(m*(p-1)+i,1) = d1/2 - e/(sqrt(2));
    else
        coor33(m*(p-1)+i,1) = (d1+d2)/4 - e/(2*sqrt(2)) + (i-(p))*e/(p*sqrt(2));
    end
    coor33(m*(p-1)+i,2) = d2/2 + e/(2*sqrt(2)) + (i-p)*e/(p*sqrt(2));
end

for i = 1:m-1 %remaining inside region
    
    for j = 1:p-1
        % p-1 because diagonal excluded 
        
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
    if mod(p,2)==0
        mult = (p+1)/2;
    else
        mult = (p+1)/2 + 1;
    end

    if i<mult - 1
        coor44(m*(p-1)+i,1) = (d1/2+e/(2*sqrt(2))) + (i)*e/(p*sqrt(2));
    elseif i == (p+2)/2 - 1
        coor44(m*(p-1)+i,1) = d1/2 + e/sqrt(2);
    else
        coor44(m*(p-1)+i,1) = (d1/2+e/(2*sqrt(2))) - (i-(p))*e/(p*sqrt(2));
    end



    %coor44(m*(p-1)+i,1) = (d1/2+e/(2*sqrt(2)))  ;
    
    
    coor44(m*(p-1)+i,2) = (d2/2-e/(2*sqrt(2))) + (i)*e/(p*sqrt(2));
    
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

%%%%   Generating elements   %%%%


EL = zeros(NoE, NPE);

for i = 1:m 
    
    for j = 1: 4*p
        
        if j == 1 %first element of the row
           
            EL((i-1)*(4*p)+j,1) = (i-1)*(4*p) + 1;
            EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) + 1;
            EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p;
            EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1;
            
        elseif j == 4*p %last element of the row
            
            EL((i-1)*(4*p)+j,1) = (i)*(4*p);%first node = last node of first row
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