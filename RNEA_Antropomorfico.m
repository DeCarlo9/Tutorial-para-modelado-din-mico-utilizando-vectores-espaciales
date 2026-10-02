%Modelado de robot de 3gdl con una fuerza externa
clc
close all
clear all

syms q1 q2 q3 q1p q2p q3p q1pp q2pp q3pp 
syms m1 m2 m3 l1 l2 l3 lc1 lc2 lc3 g I1 I2 I3 fx fy fz

f4_global = [fx; fy; fz];

%Algoritmo RNEA
%vectores de eje de giro
s1 = [0; 0; 1; 0; 0; 0];
s2 = [0; 0; 1; 0; 0; 0];
s3 = [0; 0; 1; 0; 0; 0];

%Matrices de rotación
%Matriz que referencia al marco 0 respecto al 1
R1 = [cos(q1),  sin(q1), 0; 
     -sin(q1),  cos(q1), 0; 
       0,   0, 1];

%Matriz que referencia al marco 1 respecto al 2
R2 = [cos(q2),  0, sin(q2); 
     -sin(q2),  0, cos(q2); 
            0, -1, 0];

%Matriz que referencia al marco 2 respecto al 3
R3 = [cos(q3),  sin(q3), 0; 
     -sin(q3),  cos(q3), 0; 
       0,   0, 1];

%Matriz que referencia al marco 0 respecto al 3 (misma que del marco 4 al 0)
R4 = R3* R2 * R1;

%Fuerza externa expresada en el marco 4 (efector)
f4_local = R4*f4_global;
f4 = [0;0;0;f4_local];

%vectores a centros de masa
c1 = [0; 0; lc1];
c2 = [lc2; 0; 0];
c3 = [lc3; 0; 0];

%vectores entre origenes
r1 = [0; 0; 0];
r2 = [0; 0; l1];
r3 = [l2; 0; 0];
r4 = [l3; 0; 0]; 

%Matrices de inercia local
Ic1 = [I1 0 0; 0 I1 0; 0 0 I1];
Ic2 = [I2 0 0; 0 I2 0; 0 0 I2];
Ic3 = [I3 0 0; 0 I3 0; 0 0 I3];

%condiciones iniciales
v0 = [0; 0; 0; 0; 0; 0];
a0 = [0; 0; 0; 0; 0; g];

%Matrices de transformación de movimiento
X1 = [R1,          zeros(3,3); 
     -R1 * skew(r1),   R1];

X2 = [R2,          zeros(3,3); 
     -R2 * skew(r2),   R2];

X3 = [R3,          zeros(3,3); 
     -R3 * skew(r3),   R3];

% Matriz de transformación del marco 4 al marco 3 (misma orientación)
X4 = [eye(3,3),            zeros(3,3); 
      -eye(3,3) * skew(r4), eye(3,3)];

%Recurción hacia adelante
%velocidad y aceleración espacial de eslabón 1
v1 = X1*v0 + s1*q1p;

v1x = [skew(v1(1:3)), zeros(3,3);
       skew(v1(4:6)), skew(v1(1:3))];

a1 = X1*a0 + s1*q1pp + v1x*(s1*q1p);

%velocidad y aceleración espacial de eslabón 2
v2 = X2*v1 + s2*q2p;

v2x = [skew(v2(1:3)), zeros(3,3);
       skew(v2(4:6)), skew(v2(1:3))];

a2 = X2*a1 + s2*q2pp + v2x*(s2*q2p);

%velocidad y aceleración espacial de eslabón 3
v3 = X3*v2 + s3*q3p;

v3x = [skew(v3(1:3)), zeros(3,3);
       skew(v3(4:6)), skew(v3(1:3))];

a3 = X3*a2 + s3*q3pp + v3x*(s3*q3p);

%Matrices de inercia
I1_e = [Ic1 + m1*skew(c1)*transpose(skew(c1)), m1*skew(c1);
        m1*transpose(skew(c1)), m1*eye(3,3)];

I2_e = [Ic2 + m2*skew(c2)*transpose(skew(c2)), m2*skew(c2);
        m2*transpose(skew(c2)), m2*eye(3,3)];

I3_e = [Ic3 + m3*skew(c3)*transpose(skew(c3)), m3*skew(c3);
        m3*transpose(skew(c3)), m3*eye(3,3)];

%Recurción hacia atrás
%Calculo de fuerzas espaciales
f3 = I3_e*a3 - transpose(v3x)*(I3_e*v3) - transpose(X4)*f4;
f3 = simplify(f3);

T3 = transpose(s3)*f3;
T3 = simplify(T3);

f2 = I2_e*a2 - transpose(v2x)*(I2_e*v2) + transpose(X3)*f3;
f2 = simplify(f2);

T2 = transpose(s2)*f2;
T2 = simplify(T2);

f1 = I1_e*a1 - transpose(v1x)*(I1_e*v1) + transpose(X2)*f2;
T1 = transpose(s1)*f1;
T1 = simplify(T1);
