function createGifFromTwoVectors(Y1,Y2,X,filename)
if length(X) ~= length(Y1) || length(X) ~= length(Y2)
    error('X and Y must be the same length');
end
h = figure('units','normalized','outerposition',[0 0 1 1]);
%% First subplot:
subplot(1,2,1);
box on;grid on;
xlabel('ux2','FontSize',26)
ylabel('lambda','FontSize',26)
axis tight manual;
hold on;
point1 = plot(Y1(1),X(1),'bo-','MarkerFaceColor','red');
line1  = plot(Y1(1),X(1),'b-');
axis([min(Y1),max(Y1),min(X),max(X)])
hold
%% Second subplot:
subplot(1,2,2);
box on;grid on;
xlabel('ux3','FontSize',26)
ylabel('lambda','FontSize',26)
axis tight manual;
hold on;
point2 = plot(Y2(1),X(1),'bo-','MarkerFaceColor','red');
line2  = plot(Y2(1),X(1),'b-');
axis([min(Y2),max(Y2),min(X),max(X)])
hold
%% Loop:
FR = 1:50:length(X)-10;
FR = [FR length(X)-9:length(X)];
for N = 1 : length(FR)
    n = FR(N);
    subplot(1,2,1);
    set(point1,'YData',X(n),'XData',Y1(n));
    set(line1,'YData',X(1:n),'XData',Y1(1:n));
    subplot(1,2,2);
    set(point2,'YData',X(n),'XData',Y2(n));
    set(line2,'YData',X(1:n),'XData',Y2(1:n));
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end