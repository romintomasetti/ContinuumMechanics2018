function configureFigure(fig)
figure(fig);
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultTextInterpreter','latex');

FONTSIZE  = 22;
SIZE_SAVE_X = 20;
SIZE_SAVE_Y = SIZE_SAVE_X*3/4;

grid on;hold on;
set(gca,'FontSize',FONTSIZE);
set(gcf ,'Units' , 'centimeters' )
set(gcf ,'Position' ,[0 0 SIZE_SAVE_X SIZE_SAVE_Y]);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 SIZE_SAVE_X SIZE_SAVE_Y];
box on ; grid on ; hold on ;
FONTSIZE = 22;
set(gca,'FontSize',FONTSIZE);
end