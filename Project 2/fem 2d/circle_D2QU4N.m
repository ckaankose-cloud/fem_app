function [NL1, EL1] = circle_D2QU4N(d1, d2, p, m, R, inclusion)

q = [0 0; d1 0; 0 d2; d1 d2];
q_inclusion = [d1/2-R/4 d2/2-R/4; d1/2+R/4 d2/2-R/4; d1/2-R/4 d2/2+R/4; d1/2+R/4 d2/2+R/4];

PD = 2;
if inclusion 
    NoN = 2*(p+1)*(2*m+1) + 2*(p-1)*(2*m+1) + (p-1)^2;
    NoE = 2*4*p*m+p^2;
    NPE = 4;
else
    NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1);
    NoE = 4*p*m;
    NPE = 4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Nodes %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

NL1 = zeros(NoN, PD);

a = d1/p;
b = d2/p;

a_inclusion = R/(2*p);

%%%% Region 1 - Void %%%%

if inclusion
    coor11 = zeros((p+1)*(2*m+1), PD);%15x2
else
    coor11 = zeros((p+1)*(m+1), PD); %9x2
end

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

%%%%% Region 1 - Inclusion %%%%%
if inclusion
    for i = ((p+1)*(2*m)+1):(p+1)*(2*m+1)
        coor11(i,1) = q_inclusion(1,1) + (i-1-(p+1)*(2*m))*(a_inclusion);
        coor11(i,2) = q_inclusion(1,2);
    end
    
    for i = m+1:2*m-1
        for j = 1:p+1
            dx = (coor11(2*m*(p+1)+j,1) - coor11(m*(p+1)+j,1))/m;
            dy = (coor11(2*m*(p+1)+j,2) - coor11(m*(p+1)+j,2))/m;
            
            coor11(i*(p+1)+j,1) = coor11(((i-1)*(p+1))+j,1) + dx;
            coor11(i*(p+1)+j,2) = coor11(((i-1)*(p+1))+j,2) + dy;
        end
    end
end


%%%%% Region 2 - Void  %%%%%

if inclusion
    coor22 = zeros((p+1)*(2*m+1), PD);
else
    coor22 = zeros((p+1)*(m+1), PD);
end

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

%%%%% Region 2 - Inclusion %%%%%

if inclusion
    for i = ((p+1)*(2*m)+1):(p+1)*(2*m+1)
        coor22(i,1) = q_inclusion(3,1) + (i-1-(p+1)*(2*m))*(a_inclusion);
        coor22(i,2) = q_inclusion(3,2);
    end
    
    for i = m+1:2*m-1
        for j = 1:p+1
            dx = (coor22(2*m*(p+1)+j,1) - coor22(m*(p+1)+j,1))/m;
            dy = (coor22(2*m*(p+1)+j,2) - coor22(m*(p+1)+j,2))/m;
            
            coor22(i*(p+1)+j,1) = coor22(((i-1)*(p+1))+j,1) + dx;
            coor22(i*(p+1)+j,2) = coor22(((i-1)*(p+1))+j,2) + dy;
        end
    end
end

%%%%% Region 3 - Void %%%%%

if inclusion
    coor33 = zeros((p-1)*(2*m+1), PD);
else
    coor33 = zeros((p-1)*(m+1), PD);
end

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

%%%%% Region 3 - Inclusion %%%%%

if inclusion
    for i = ((p-1)*(2*m)+1):(p-1)*(2*m+1)
        coor33(i,1) = q_inclusion(1,1);
        coor33(i,2) = q_inclusion(1,2) + (i-(p-1)*(2*m))*(a_inclusion);
    end
    
    for i = m+1:2*m-1
        for j = 1:p-1
            dx = (coor33(2*m*(p-1)+j,1) - coor33(m*(p-1)+j,1))/m;
            dy = (coor33(2*m*(p-1)+j,2) - coor33(m*(p-1)+j,2))/m;
            
            coor33(i*(p-1)+j,1) = coor33(((i-1)*(p-1))+j,1) + dx;
            coor33(i*(p-1)+j,2) = coor33(((i-1)*(p-1))+j,2) + dy;
        end
    end
end

%%%%% Region 4 - Void %%%%%

if inclusion
    coor44 = zeros((p-1)*(2*m+1), PD);
else
    coor44 = zeros((p-1)*(m+1), PD);
end

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

%%%%% Region 4 - Inclusion %%%%%

if inclusion
    for i = ((p-1)*(2*m)+1):(p-1)*(2*m+1)
        coor44(i,1) = q_inclusion(2,1);
        coor44(i,2) = q_inclusion(2,2) + (i-(p-1)*(2*m))*(a_inclusion);
    end
    
    for i = m+1:2*m-1
        for j = 1:p-1
            dx = (coor44(2*m*(p-1)+j,1) - coor44(m*(p-1)+j,1))/m;
            dy = (coor44(2*m*(p-1)+j,2) - coor44(m*(p-1)+j,2))/m;
            
            coor44(i*(p-1)+j,1) = coor44(((i-1)*(p-1))+j,1) + dx;
            coor44(i*(p-1)+j,2) = coor44(((i-1)*(p-1))+j,2) + dy;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if inclusion 
    for i = 1:(1+inclusion)*m+1
        NL1((i-1)*4*p+1:(i)*4*p,:) = [coor11((i-1)*(p+1)+1:(i)*(p+1),:);
        coor44((i-1)*(p-1)+1:(i)*(p-1),:);
        flipud(coor22((i-1)*(p+1)+1:(i)*(p+1),:));
        flipud(coor33((i-1)*(p-1)+1:(i)*(p-1),:))];
    end
else
    for i = 1:m+1
        NL1((i-1)*4*p+1:(i)*4*p,:) = [coor11((i-1)*(p+1)+1:(i)*(p+1),:);
        coor44((i-1)*(p-1)+1:(i)*(p-1),:);
        flipud(coor22((i-1)*(p+1)+1:(i)*(p+1),:));
        flipud(coor33((i-1)*(p-1)+1:(i)*(p-1),:))];
    end
end

%%% adding inner nodes into the NL %%%
if inclusion 
    a = R/(2*p);
    n=2*(p+1)*(2*m+1) + 2*(p-1)*(2*m+1);
    for i = 1:p-1
        for j = 1:p-1
            NL1(n+1,1) = q_inclusion(1,1) + (j)*a;
            NL1(n+1,2) = q_inclusion(1,2) + (i)*a;
            n = n+1;
        end   
    end
end


%%%%%%%%%%%%%%%%%%%%
%%%%% Elements %%%%%
%%%%%%%%%%%%%%%%%%%%

EL1 = zeros(NoE, NPE);

for i = 1:(inclusion + 1)*m
    for j = 1:4*p
        if j == 1
            EL1((i-1)*(4*p)+j,1) = (i-1)*(4*p) + 1;
            EL1((i-1)*(4*p)+j,2) = EL1((i-1)*(4*p)+j,1) + 1;
            EL1((i-1)*(4*p)+j,4) = EL1((i-1)*(4*p)+j,1) + 4*p;
            EL1((i-1)*(4*p)+j,3) = EL1((i-1)*(4*p)+j,4) + 1;
        elseif j == 4*p
            EL1((i-1)*(4*p)+j,1) = (i)*(4*p);
            EL1((i-1)*(4*p)+j,2) = (i-1)*(4*p) + 1;
            EL1((i-1)*(4*p)+j,3) = EL1((i-1)*(4*p)+j,1) + 1;
            EL1((i-1)*(4*p)+j,4) = EL1((i-1)*(4*p)+j,1) + 4*p;
        else
            EL1((i-1)*(4*p)+j,1) = EL1((i-1)*(4*p)+j-1,2);
            EL1((i-1)*(4*p)+j,4) = EL1((i-1)*(4*p)+j-1,3);
            EL1((i-1)*(4*p)+j,3) = EL1((i-1)*(4*p)+j,4) + 1;
            EL1((i-1)*(4*p)+j,2) = EL1((i-1)*(4*p)+j,1) + 1;
        end
    end
end

if inclusion
    El_sqr = zeros((p)^2,NPE);
    
        for j=1:4*(p-1)
            
            if j<p
                El_sqr(j,1) = j;
                El_sqr(j,2) = j+1;
                El_sqr(j,3) = 4*p+j;
                El_sqr(j,4) = 4*p+j-1;
            elseif j==p
                El_sqr(j,1) = j;
                El_sqr(j,2) = j+1;
                El_sqr(j,3) = j+2;
                El_sqr(j,4) = 4*p+j-1;
            elseif j>=p && j<2*p-1
                El_sqr(j,1) = j+1;
                El_sqr(j,2) = j+2;
                El_sqr(j,3) = 4*p+(p-1)*(j-p+1);
                El_sqr(j,4) = 4*p+(p-1)*(j-p);
            elseif j==2*p-1
                El_sqr(j,2) = j+1;
                El_sqr(j,3) = j+2;
                El_sqr(j,4) = j+3;
                El_sqr(j,1) = (p+1)^2;
            elseif j>2*p-1 && j<3*(p-1)+1
                El_sqr(j,1) = j+2;
                El_sqr(j,2) = j+3;
                El_sqr(j,3) = (p+1)^2-(j-2*p)-1;
                El_sqr(j,4) = (p+1)^2-(j-2*p); 
            elseif j==3*(p-1)+1
                El_sqr(j,3) = j+2;
                El_sqr(j,4) = j+3;
                El_sqr(j,1) = j+4;
                El_sqr(j,2) = (p+1)^2-(p-2);
            else
                El_sqr(j,1) = j+3;
                El_sqr(j,2) = j+4;
                El_sqr(j,3) = (p+1)^2-(j-3*(p-1))*(p-1)+1;
                El_sqr(j,4) = (p+1)^2-(j-3*(p-1)-1)*(p-1)+1;
                
            end
        end
        
        for i=1:p-2
            for j=1:p-2
                El_sqr(4*(p-1)+(i-1)*(p-2)+j,1) = 4*p+(i-1)*(p-1)+j;
                El_sqr(4*(p-1)+(i-1)*(p-2)+j,2) = 4*p+(i-1)*(p-1)+j+1;
                El_sqr(4*(p-1)+(i-1)*(p-2)+j,3) = 4*p+(i-1)*(p-1)+j+p;
                El_sqr(4*(p-1)+(i-1)*(p-2)+j,4) = 4*p+(i-1)*(p-1)+j+p-1;
            end
        end
    
    El_sqr = El_sqr + size(EL1,1) - p^2;
    EL1(NoE-p^2+1:NoE,:) = El_sqr(:,:);
end

end
