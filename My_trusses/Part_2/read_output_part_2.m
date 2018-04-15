%function read_output_part_2
close all
direct = ['D:\GOOGLE DRIVE\UNIVERSITY\MASTER 1\SECOND_QUADRIMESTER\Continuum Mechanics'...
    '\Project\My_trusses\Part_2\'];

ancienne_dir = pwd;
cd(direct);

%% SPHERICAL_ARC_LENGTH:
files = {'SPHERICAL_ARC_LENGTH\Node_2_DISPLACEMENTS.ascii',...
    'SPHERICAL_ARC_LENGTH\Node_3_DISPLACEMENTS.ascii',...
    'SPHERICAL_ARC_LENGTH\Bar_2_internal_force.ascii',...
    'SPHERICAL_ARC_LENGTH\Lambda.ascii'};
vector_names = {'ux2','uy2','ux3','uy3',...
                'Fint1','Fint2','Fint3','Fint4',...
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
plot_struct(spherical)
%% 
%end