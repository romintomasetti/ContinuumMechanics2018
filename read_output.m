function read_output()
close all
direct = ['D:\GOOGLE DRIVE\UNIVERSITY\MASTER 1\SECOND_QUADRIMESTER\Continuum Mechanics'...
    '\Project\My_trusses\tests\'...
    'test2bar_updated_normal_plane_arc_length\'];

name = 'updated_normal_plane_arc_length';

ancienne_dir = pwd;
cd(direct);

%% Parameters:
a = 2.0;
b = 1.0;
l0 = sqrt(a^2 + b^2);
E = 70e9;
A = 0.01;

counter = 1;
lower_bound = 0;
upper_bound = 2.2;
step        = 0.001;
u_b = zeros(1,(upper_bound-lower_bound)/step);
for n = lower_bound : step : upper_bound
    u_b(counter) = n;
    counter = counter + 1;
end
P = (E*A*b*b*b)/(2*l0*l0*l0) .* u_b .* (u_b-1) .* (u_b-2);
u1 = (1+sqrt(3)/3)*b;
u2 = (1-sqrt(3)/3)*b;
Pcrit_analytical = E*A/(2*l0^3)*(u2^3-3*b*u2^2+2*b^2*u2);
%%
f = fopen('Node_2_DISPLACEMENTS.ascii', 'r');
uy2 = [];
while ~feof(f)
    line = fgetl(f);
    R = sscanf(line,'%f %f');
    uy2 = [uy2 R(2)];
end
fclose(f);
%%
f = fopen('Lambda.ascii','r');
lambda_ = [];
while ~feof(f)
    line = fgetl(f);
    R = sscanf(line,'%f');
    lambda_ = [lambda_ R(1)];
end
fclose(f);

%%
X = -uy2/b;
Y = 1.5*Pcrit_analytical*lambda_;

size(X)
size(Y)

cd(ancienne_dir);
configureFigure(figure);
cd(direct);
hold on;
plot(X,Y(1:end),'bs','linewidth',3.0);
plot(u_b,P,'r-.','linewidth',3.0);
xlabel('p/b[-]')
ylabel('q[N]')
leg = legend(['Steps : ' num2str(length(X))]);
leg.Location = 'best';
saveas(gcf,['AnalyticalPART1_' name '.eps'],'epsc2');

%% Compute the error:
poly = polyfit(X,Y,10);
for i = 1 : length(P)
   abs_error(i) = abs(polyval(poly,u_b(i))-P(i));
end
figure;
plot(u_b,abs_error,'o-')
int_of_error = cumsum((abs_error(2:end)+abs_error(1:end-1))/...
                    2*step);
                
txt = sprintf('%0.10e',int_of_error(end))
fprintf('trapz of the abs. error: %s.\n', txt);

%%
cd(ancienne_dir);
return