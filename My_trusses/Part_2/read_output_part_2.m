%function read_output_part_2
close all
direct = ['D:\GOOGLE DRIVE\UNIVERSITY\MASTER 1\SECOND_QUADRIMESTER\Continuum Mechanics'...
    '\ContinuumMechanics2018\My_trusses\Part_2\'];

ancienne_dir = pwd;
cd(direct);

%% SPHERICAL_ARC_LENGTH:
method = 'UPDATED_NORMAL_PLANE';
method = 'SPHERICAL_ARC_LENGTH';
method = 'NEWTON_RAPHSON';
files = {[method '\Node_2_DISPLACEMENTS.ascii'],...
    [method '\Node_3_DISPLACEMENTS.ascii'],...
    [method '\Lambda.ascii']};
vector_names = {'ux2','uy2','ux3','uy3',...
                'lambda'};
spherical.INFO = 'Contains everything for spherical a.-l.';
counter_vector = 1;
for i = 1 : length(files)
    f = fopen(files{i}, 'r');
    line = fgetl(f);
    R = sscanf(line,'%f %f');
    %disp([files{i} ' size ' num2str(size(R))])
    for j = 1 : size(R,1)
        spherical.(vector_names{counter_vector}) = [];
        counter_vector = counter_vector + 1;
    end
    frewind(f);
    while ~feof(f)
        line = fgetl(f);
        R = sscanf(line,'%f %f');
        for j = 1 : size(R,1)
            spherical.(vector_names{counter_vector-j}) = ...
                [spherical.(vector_names{counter_vector-j*1}) R(length(R)-j+1)];
        end
    end
    fclose(f);
end
spherical
comet(spherical.uy2,spherical.lambda)
%plot_struct(spherical)
% createGifFromTwoVectors(spherical.ux2,...
%     spherical.uy2,spherical.lambda,'test.gif')
% configureFigure(figure)
% plot(spherical.lambda,'*-')
% xlabel('step')
% ylabel('lambda')
%% 
%end