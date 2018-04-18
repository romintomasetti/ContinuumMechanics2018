function plot_struct(given,method,save)
disp('>> plot_struct(begin)');
%% Get the field names:
names = fieldnames(given)
given
nbr   = length(names);

%% Plot fields 2 and 3 together:
configureFigure(figure);
plot(given.(names{2}),given.lambda,'ro','MarkerFaceColor','red',...
    'MarkerSize',2)
xlabel('$u_{(2,x)}$')
ylabel('$\lambda$')
if save == 1
    saveas(gcf,[method '_ux2.eps'],'epsc2');
end

configureFigure(figure);
plot(given.(names{3}),given.lambda,'ro','MarkerFaceColor','red',...
    'MarkerSize',2)
xlabel('$u_{(2,y)}$')
ylabel('$\lambda$')
if save == 1
    saveas(gcf,[method '_uy2.eps'],'epsc2');
end

%% Plot fields 4 and 5 together:
% configureFigure(figure);
% plot(given.(names{4}),given.lambda,'ro','MarkerFaceColor','red',...
%     'MarkerSize',2)
% xlabel('$u_{(3,x)}$')
% ylabel('$\lambda$')
% saveas(gcf,[method '_ux3.eps'],'epsc2')
% 
% configureFigure(figure);
% plot(given.(names{5}),given.lambda,'ro','MarkerFaceColor','red',...
%     'MarkerSize',2)
% xlabel('$u_{(3,y)}$')
% ylabel('$\lambda$')
% saveas(gcf,[method '_uy3.eps'],'epsc2')
%% Plot internal force:
% ampl1 = sqrt((given.Fint1).^2+(given.Fint2).^2);
% ampl2 = sqrt((given.Fint3).^2+(given.Fint4).^2);
% figure;
% subplot(121)
% hold on;
% plot(given.Fint1);
% plot(given.Fint2);
% legend('Fint1','Fint2');
% subplot(122);
% hold on;
% plot(given.Fint3);
% plot(given.Fint4);
% legend('Fint3','Fint4');

%% Animation:
%return
a     = 0.75
b     = 0.25
alpha = 30*pi/180


figure;

node1 = [0,0];
node4 = [2*a*cos(alpha)-b,0];

node2 = [a*cos(alpha)  , a*sin(alpha)];
node3 = [a*cos(alpha)-b, a*sin(alpha)];

xData = [node1(1),node2(1),node3(1),node4(1)];
yData = [node1(2),node2(2),node3(2),node4(2)];

axis([node1(1),node4(1),-a,a]);
hold;
bars_plot  = plot(xData,yData,'o-','Color','black','MarkerFaceColor','red',...
    'linewidth',3);

for frame = 1 : length(given.(names{2}))
    % Update bar1:
    xData(2) = node2(1) + given.(names{2})(frame);
    yData(2) = node2(2) + given.(names{3})(frame);
    % Update bar2:
    xData(3) = node3(1) + given.(names{4})(frame);
    yData(3) = node3(2) + given.(names{5})(frame);
    % Update bar3:
    set(bars_plot,'XData',xData,'YData',yData);
    %title(['Step ' num2str(frame) ' over ' num2str(length(given.(names{2})))]);
    drawnow;
end

end