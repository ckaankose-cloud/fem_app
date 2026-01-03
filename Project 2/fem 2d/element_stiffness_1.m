function k= element_stiffness_1(nl, ENL, E, A) %takes input as nodelist nl (connectivity per element) and enl 

X(1)= ENL(nl(1),1);%ENL(1,1) first node x 
Y(1)= ENL(nl(1),2);%ENL(1,2) first node y
X(2)= ENL(nl(2),1);%ENL(2,1) second node x   
Y(2)= ENL(nl(2),2);%ENL(2,2) second node y
    
L=sqrt((X(1)-X(2))^2 + (Y(1)-Y(2))^2); %length 

C= (X(2)-X(1))/L;
S= (Y(2)-Y(1))/L;

k= (A*E/L)* [C^2 C*S -C^2 -C*S
             C*S S^2 -C*S -S^2
             -C^2 -C*S C^2 C*S
             -C*S -S^2 C*S S^2];
end