clear all
close all
clc;

%NPE : Number of Nodes per Element
%%%%%% PRE-PROCESS %%%%%%

format long

d1 = 1;
d2 = 1;
p = 3;
m = 4;
R = 0.2;
element_type = 'D2TR4N'; % 'D2TR3N'

[NL, EL] = void_mesh(d1, d2, p, m, R, element_type);

%%%%%%%%%%%%% PROCESS %%%%%%%%%%%%%%%

BC_type = 'Extension'; % Expansion-Shear

[ENL, DOFs, DOCs] = assign_BCs(NL, BC_type);

K = assemble_stiffness(ENL, EL, NL);

Fp = assemble_forces(ENL, NL);

Up = assemble_displacements(ENL, NL);

Kpu = K(1:DOFs,1:DOFs);
Kpp = K(1:DOFs,DOFs+1:DOFs+DOCs);
Kuu = K(DOFs+1:DOCs+DOFs,1:DOFs);
Kup = K(DOFs+1:DOCs+DOFs,DOFs+1:DOCs+DOFs);

F = Fp - Kpp * Up;

Uu = inv(Kpu)*F;

Fu = Kuu*Uu + Kup*Up;

ENL = update_nodes(ENL,Uu,NL,Fu);

post_process(NL, EL, ENL)
