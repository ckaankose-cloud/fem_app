clear all;
close all; 
clc;

format long; 
d1 = 1;
d2= 1;
p =3;
m = 3;
R = 0.2;
element_type = "D2QU4N"; %"D2TR3N"

[NL, EL] = void_mesh(d1, d2, p, m, R, element_type);

%PROCESS%
[ENL,DOFs,DOCs]= assign_BCs(NL);%takes NL input and outputs ENL,DOFs,DOCs

K= assemble_stiffness(ENL,EL,NL);
Fp=assemble_forces(ENL, NL);
Up=assemble_displacements(ENL,NL);  
%assemble forces and displacements 
Kpu= K(1:DOFs,1:DOFs); %DOFxDOF                                                                         m kü
Kpp= K(1:DOFs, DOFs+1: DOFs+DOCs); %DOFxDOC
Kuu= K(DOFs+1:DOCs+DOFs,1:DOFs);% DOCxDOF
Kup= K(DOFs+1:DOCs+DOFs,DOFs+1:DOCs+DOFs);%DOCxDOC

F = Fp - Kpp * Up;

Uu=inv(Kpu)*F;  

Fu=Kuu*Uu+Kup*Up;

ENL= update_nodes(ENL,Uu,NL,Fu);

Node_flag= 'on';    %shows node number 
Element_flag= 'on'; %shows element number 
mag=1;%maginification factor exaggarates the small deformation  
post_process(NL,EL,ENL)