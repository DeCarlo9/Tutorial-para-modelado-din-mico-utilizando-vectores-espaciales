%CRBA del robot antropomórfico
clc
close all
clear all

syms m1 m2 m3 real
syms lc1 lc2 lc3 real
syms l1 l2 l3 real
syms I1 I2 I3 real
syms q1 q2 q3 real

% Ejes de rotación 
s1 = [0; 0; 1; 0; 0; 0];
s2 = [0; 0; 1; 0; 0; 0];
s3 = [0; 0; 1; 0; 0; 0];

% Transformaciones locales espaciales
X1 = rotz(q1);
X2 = rotz(q2) * rotx(sym(pi)/2) * xlt([0; 0; l1]); 
X3 = rotz(q3) * xlt([l2; 0; 0]);
X4 = xlt([l3; 0; 0]); 

% Matrices de Inercia Espacial
I1_e = mcI(m1, [0; 0; lc1], diag([I1, I1, I1]));
I2_e = mcI(m2, [lc2; 0; 0], diag([I2, I2, I2]));
I3_e = mcI(m3, [lc3; 0; 0], diag([I3, I3, I3]));

% Etapa 1: Inercias del cuerpo compuesto (propogación hacia atrás)
Ic3 = I3_e;
Ic2 = I2_e + X3' * Ic3 * X3;
Ic1 = I1_e + X2' * Ic2 * X2;

% Cálculo de fuerzas espaciales compuestas
F3 = Ic3 * s3;
F2 = Ic2 * s2;
F1 = Ic1 * s1;

% Etapa 2: Construcción de la matriz de inercia 
M = sym(zeros(3,3));

M(3,3) = s3' * F3;

F23 = (X3' * F3);

M(2,3) = s2' * F23;
M(3,2) = M(2,3);

F13 = (X2' * X3' * F3);
M(1,3) = s1' * F13;
M(3,1) = M(1,3);

M(2,2) = s2' * F2;

F12 = (X2' * F2);
M(1,2) = s1' * F12; 
M(2,1) = M(1,2);

M(1,1) = s1' * F1;

M = simplify(M);