clear all
close all
clc;

%%%%%% PRE-PROCESS %%%%%%

format long

d1 = 1;

d2 = 1;

p = 4;

m = 4;

R = 0.2;

element_type = 'D2QU4N';  %'D2TR3N' 'D2TR6N' 'D2QU4N' 'D2QU8N' 'D2QU9N'

inc = 0;

PR = [1 0.3; 10 0.3];

magnitude = 0.1;

[NL, EL] = void_mesh(d1, d2, p, m, R, element_type, inc);
%[NL, EL] = uniform_mesh(d1, d2, p, m, element_type);

%%%%%%%%%%%%% PROCESS %%%%%%%%%%%%%%%

BC_type = 'Extension'; % Extension, Shear, Expansion

[ENL, DOFs, DOCs] = assign_BCs(NL, BC_type, magnitude, d1, d2);

K = assemble_stiffness(ENL, EL, NL, PR, p, m, element_type);

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

post_process(NL, EL, ENL,PR,p,m,element_type);


