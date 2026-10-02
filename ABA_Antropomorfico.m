%ABA del robot antropomórfico
clc
clear all

syms q1 q2 q3 q1p q2p q3p real
syms T1 T2 T3 m1 m2 m3 I1 I2 I3 real
syms l1 l2 l3 lc1 lc2 lc3 g fx fy fz real

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

% Condición inicial                                          
v0 = zeros(6,1);
a0 = [0; 0; 0; 0; 0; g]; 

% Matrices de rotación
R1 = X1(1:3, 1:3);
R2 = X2(1:3, 1:3);
R3 = X3(1:3, 1:3);

% Matriz de rotación del marco 4 al marco 0
R_0_4 = R3 * R2 * R1; 
f_loc4 = R_0_4 * [fx; fy; fz];
f4 = [0; 0; 0; f_loc4];

% Pasada 1: Cinemática hacia adelante
v1 = X1*v0 + s1*q1p;
c1 = crm(v1)*(s1*q1p);
p1 = crf(v1)*(I1_e*v1);

v2 = X2*v1 + s2*q2p;
c2 = crm(v2)*(s2*q2p);
p2 = crf(v2)*(I2_e*v2);

v3 = X3*v2 + s3*q3p;
c3 = crm(v3)*(s3*q3p);
p3 = crf(v3)*(I3_e*v3);

% Pasada 2: Inercia articulada hacia atrás
I3_A = I3_e;
p3_A = p3 - X4'*f4; 
U3 = I3_A*s3;
D3 = s3'*U3;
u_3 = T3 - s3'*p3_A;
I3a = I3_A - (U3*U3')/D3;
p3a = p3_A + I3a*c3 + U3*(u_3/D3);

I2_A = I2_e + X3'*I3a*X3;
p2_A = p2 + X3'*p3a;
U2 = I2_A*s2;
D2 = s2'*U2;
u_2 = T2 - s2'*p2_A;
I2a = I2_A - (U2*U2')/D2;
p2a = p2_A + I2a*c2 + U2*(u_2/D2);

I1_A = I1_e + X2'*I2a*X2;
p1_A = p1 + X2'*p2a; 
U1 = I1_A*s1;
D1 = s1'*U1;
u_1 = T1 - s1'*p1_A;

% Pasada 3: Aceleraciones articuladas hacia adelante
a1_p = X1*a0 + c1;
q1pp = (u_1 - U1'*a1_p)/D1;
a1   = a1_p + s1*q1pp;

a2_p = X2*a1 + c2;
q2pp = (u_2 - U2'*a2_p)/D2;
a2   = a2_p + s2*q2pp;

a3_p = X3*a2 + c3;
q3pp = (u_3 - U3'*a3_p)/D3;
