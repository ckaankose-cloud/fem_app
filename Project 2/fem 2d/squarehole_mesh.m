function [NL, EL, NoE, NoN] = squarehole_mesh(p_input)

d1 = 1;
d2 = 1;
q = [0 0; d1 0; 0 d2; d1 d2];
m = 4;
R = 0;
p=p_input;

PD = 2;
NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1);

%%%%
NoE = 16*p*m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Nodes   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NL = zeros(NoN, PD);

a = d1/p;
b = d2/p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Region 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Region 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Region 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Region 4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%karnin nodeları
for i = 1 : p*m
    NL(NoN,:)=[];
    NoN = NoN -1;
end


%buraya kadar nodeları sıraladı
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Elements   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EL = [];
%her bir kareyi çiziyor(ilk aşama)

a=[];
b=[];
for i = 1: m
    for k = (i-1)*p*m+1 : i*p*m-1
        a = [k,k+1];
        b =[b;a];
        
    end
    b = [b; i*p*m i*p*m-p*m+1];
end
%b

%her bir diagonalı çiziyor

c =[];
d= [];
NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1);
for i = 1 : NoN-p*m-p*m
    c = [i, i+ p*m];
    d =[d;c];
end
%d


e=[];
f=[];
for i = 1: m
    for k = (i-1)*p*m+1 : i*p*m-1
        e = [k,k+p*m+1];
        f =[f;e];
        
    end
    f = [f; i*p*m i*p*m+1];
end
%f
EL=[EL;b];
EL=[EL;d];
EL=[EL;f];
for i = 1 : p*m
    EL(end,:)=[];    
end
NoN = size(NL,1);
NoE = size(EL,1);
