function read_output_newton_in_time
%% Open files:
f_infos = fopen('Some_more_infos.ascii');
f_node2 = fopen('Node_2_DISPLACEMENTS.ascii');
f_PK2   = fopen('NR_in_time_PK2.ascii');

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

%% Read the PK2 stress inside the bar:
counter = 1;
while ~feof(f_PK2)
    line = fgetl(f_PK2);
    R = sscanf(line,'t=%f;PK2_bar%d=%f');
    PK2(counter) = R(3);
    counter = counter + 1;
end

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
%comet(uy2,F);
plot(uy2,F);
xlabel('uy2')
ylabel('loading')

%% Plot the structure in time:
figure;
subplot(1,2,1);
box on; grid on;
axis([0,max(t),min(F),max(F)]);
hold;
time    = 0;
loading = 0;
h2 = plot(time,loading,'.-','MarkerSize',5);
subplot(1,2,2);
box on; grid on;
a = 2;
b = 1;
X = [0,a,2*a];
Y = [0,b,0];

axis([0,2*a,b+min(uy2),b+max(uy2)]);
hold;

h1 = plot(X,Y,'o-');
for i = 1 : length(t)
    Y(2) = b + uy2(i);
    set(h1,'YData',Y);
    time    = [time    t(i)];
    loading = [loading F(i)];
    set(h2,'XData',time,'YData',loading);
    drawnow;
    pause(0.05);
end

%% Plot the stresses:
figure;
subplot(1,2,1);
plot(t,PK2,'.-');
xlabel('time [???]');
ylabel('PK2 stress [???]');

subplot(1,2,2);
plot(uy2,PK2,'.-');
xlabel('uy2 [??]');
ylabel('PK2 stress [???]');


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