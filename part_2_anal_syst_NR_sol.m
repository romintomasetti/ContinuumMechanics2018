%% Newton-Raphson solution for the analytical system of truss 2

syms u v q
syms eq1(u,v,q)
syms eq2(u,v,q)
syms deriv1u(u,v,q)
syms deriv1v(u,v,q)
syms deriv2u(u,v,q)
syms deriv2v(u,v,q)
syms Jacobian(u,v,q)
syms Inv_Jac(u,v,q)

q_max = 23962861.8909;
E1 = 70e9;
A1 = 25e-4;
a = 0.75;
b = 0.25;
theta = 30*pi/180;
E2 = E1;
A2 = A1;


eq1(u,v,q) = E1*A1/a^3*(2*a*cos(theta)*u+u^2-2*v*a*sin(theta)+v^2)*...
    (a*cos(theta)+u) + A2*E2/b^3*(b*u+u^2)*(b/2+u)*8;

eq2(u,v,q) = E1*A1/a^3*(2*cos(theta)*a*u+u^2-...
    2*sin(theta)*a*v+v^2)*(v-a*sin(theta)) - 2 *q
    
deriv1u(u,v,q) = diff(eq1, u);
deriv1v(u,v,q) = diff(eq1, v);
deriv2u(u,v,q) = diff(eq2, u);
deriv2v(u,v,q) = diff(eq2, v);
disp('Derivative computed !')

Jacobian(u,v,q) = [deriv1u(u,v,q) deriv1v(u,v,q);
                deriv2u(u,v,q) deriv2v(u,v,q)];
            
Inv_Jac(u,v,q)  = inv(Jacobian)

disp('Inverse Jacobian computed !')
            
startingPoint = [0;0];

lambda = 1e-3;
dlambda= 1e-3;

XY = [0,0];
lambda_vec = [lambda];

configureFigure(figure);
subplot(1,2,1);
hold on;
pl = plot(XY(:,1),lambda_vec,'r-');
%plot(spherical.ux2(1:100),spherical.lambda(1:100));
xlabel('$u$');
ylabel('$\lambda$');
subplot(1,2,2)
pl2 = plot(XY(:,2),lambda_vec,'r-');
xlabel('$v$');
ylabel('$\lambda$');

for i=1:1000
    q_curr = lambda * q_max;
    inv_jac_num = Inv_Jac(startingPoint(1),startingPoint(2),q_curr);
    inv_jac_num = double(inv_jac_num);
    newPoint = startingPoint - ...
        inv_jac_num*[double(eq1(startingPoint(1),startingPoint(2),q_curr));
                     double(eq2(startingPoint(1),startingPoint(2),q_curr))];
    % Evaluate error:
    eq1_num = double(eq1(newPoint(1),newPoint(2),q_curr));
    eq2_num = double(eq2(newPoint(1),newPoint(2),q_curr));
    error = sqrt(eq1_num^2+eq2_num^2);
    iter = 0;
    while error > 1e-6
        % Continue to iterate:
        inv_jac_num = Inv_Jac(startingPoint(1),startingPoint(2),q_curr);
        inv_jac_num = double(inv_jac_num);
        newPoint = startingPoint - ...
            inv_jac_num*[double(eq1(startingPoint(1),startingPoint(2),q_curr));
                         double(eq2(startingPoint(1),startingPoint(2),q_curr))];
        % Evaluate error:
        eq1_num = double(eq1(newPoint(1),newPoint(2),q_curr));
        eq2_num = double(eq2(newPoint(1),newPoint(2),q_curr));
        error = sqrt(eq1_num^2+eq2_num^2);
        fprintf('\t>> %f\n',error);
        iter = iter + 1;
        startingPoint = newPoint;
    end
    fprintf('>> Error at step %d is %f.\n',i,error);
    lambda_vec(i) = lambda;
    
    startingPoint = newPoint;
    XY(i,:) = newPoint;
    if mod(i,1) == 0
        set(pl,'XData',XY(:,1),'YData',lambda_vec)
        set(pl2,'XData',XY(:,2),'YData',lambda_vec)        
    end
    drawnow;
    lambda = lambda + dlambda;
end
save('sol_anal_sys_NR','XY')