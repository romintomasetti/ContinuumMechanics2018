q_max = 23962861.8909;
E1 = 70e9;
A1 = 25e-4;
a = 0.75;
b = 0.25;
theta = 30*pi/180;
E2 = E1;
A2 = A1;

U_python = spherical.ux2;
V_python = spherical.uy2;

L_python = spherical.lambda;

eq1 = @(u,v,q) 2*q + E1*A1 * (u^2 + v^2 +...
        2*a*(u*cos(theta) - v*sin(theta)) *...
        (a*sin(theta) - v))/(a^3);
    
eq2 = @(u,v,q) (16*E2*A2*u*(b+u)*(b+2*u))/(b^3) +...
        E1*A1 * (u^2 + v^2 +...
        2*a*(u*cos(theta) - v*sin(theta)) *...
        (a*cos(theta) + u))/(a^3);
    
eq1(U_python(1),V_python(1),L_python(1)*q_max)
eq2(U_python(1),V_python(1),L_python(1)*q_max)