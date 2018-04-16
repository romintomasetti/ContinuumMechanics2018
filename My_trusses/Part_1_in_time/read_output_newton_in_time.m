function read_output_newton_in_time
%% Open files:
f_infos = fopen('Some_more_infos.ascii');
f_node2 = fopen('Node_2_DISPLACEMENTS.ascii');

%% Read info file:
counter = 1;
while ~feof(f_infos)
    line = fgetl(f_infos);
    R = sscanf(line,'t=%f;dt=%f;F_max=%f;T=%f');
    t(counter) = R(1);
    if(counter == 1)
        dt    = R(2);
        F_max = R(3);
        T     = R(4);
    end
    counter = counter + 1;
end
fclose(f_infos);

%% Read displacement of node 2:
counter = 1;
while ~feof(f_node2)
    line = fgetl(f_node2);
    R = sscanf(line,'%f %f');
    uy2(counter) = R(2);
    counter = counter + 1;
end
fclose(f_node2);

%% Compute the applied load:
F = zeros(1,length(t));
for i = 1 : length(F)
    F(i) = get_sawtooth(t(i),F_max,T);
end

%% Plot the applied load and the incurred displacement:
figure;

subplot(1,3,1);
plot(t,F,'o-');
xlabel('time');
ylabel('loading [???]');

subplot(1,3,2);
plot(t,uy2,'o-');
xlabel('time');
ylabel('uy2 [???]');

subplot(1,3,3);
comet(uy2,F)
xlabel('uy2')
ylabel('loading')

end

function F = get_sawtooth(t,F_max,T)
t_equ = t/T - floor(t/T);
if(t_equ <= 0.25)
    F = 4.0*F_max*t_equ;
elseif (t_equ > 0.25 && t_equ <= 0.75)
    F = 2.0*F_max*(-2.0*t_equ+1.0);
else
    F = 4.0*F_max*(t_equ-1.0);
end
end