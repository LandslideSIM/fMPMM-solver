clear,close

GIMP         = load('.\GIMP_solution\data\GIMP_data_1good.mat');
data_surface = load('Data_Buietal_2008_experimental_surface.txt');
data_failure = load('Data_Buietal_2008_experimental_failure_surface.txt');

x = GIMP.mpD.x(:,1)+0.4;
y = GIMP.mpD.x(:,2)    ;
epII = GIMP.mpD.epII;
epII(epII>0.01)=1
epII(epII<0.01)=0

fig1=figure(1)
set(fig1,'Units','pixels','Position',[100 100 541 429]);
hold on
ax3=scatter(x,y,10,epII,'filled');
map = [0 1 0;1 0 0]; colormap(map);
pp.cbchar='Region';
pp.cbpos =[0.32 0.75 0.4 0.05];
pos            = pp.cbpos;
pp.cblpos=[pos(1) pos(2)+2];

cb=colorbar('FontSize',12,'TickLabelInterpreter','latex','Fontsize',12,'Location','SouthOutside');
cb.Position         =pp.cbpos;
cb.Label.String     =pp.cbchar;
cb.Label.FontSize   =12;
cb.Label.Interpreter='Latex';
%cb.Label.Position   =pp.cblpos;
cb.Ticks      = [0.25 0.75];
cb.TickLabels = {'Undeformed','Deformed'};

ax1=plot(data_surface(:,1),data_surface(:,2),'bs:','LineWidth',3);
ax2=plot(data_failure(:,1),data_failure(:,2),'bo:','LineWidth',3);
box on;
tit = {'Experiment: final geometry','Experiment: failure surface','Numerical solution'};
h1=legend([ax1 ax2 ax3],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.4451 0.5412 0.4361 0.1490],'NumColumns',1);
hold off;
axis equal;
xlim([0 0.5]);
ylim([0 0.2]);
xlabel('$x$ (mm)');
ylabel('$y$ (mm)');
grid on;
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
% print(gcf,'figElastoPlasticCollapse','-depsc');
% print(gcf,'figElastoPlasticCollapse','-dpng');
% 
% 

fig2=figure(2)
set(fig2,'Units','pixels','Position',[100 100 541 429]);
epII = GIMP.mpD.epII;
scatter(x,y,10,epII,'filled');
colormap(jet);
caxis([0 8])
pp.cbchar='Region';
pp.cbpos =[0.32 0.75 0.4 0.05];
pos            = pp.cbpos;
pp.cblpos=[pos(1) pos(2)+2];

cb=colorbar('FontSize',12,'TickLabelInterpreter','latex','Fontsize',12,'Location','SouthOutside');
cb.Position         =pp.cbpos;
cb.Label.String     =pp.cbchar;
cb.Label.FontSize   =12;
cb.Label.Interpreter='Latex';
%cb.Label.Position   =pp.cblpos;
box on;
axis equal;
xlim([0 0.5]);
ylim([0 0.2]);
xlabel('$x$ (mm)');
ylabel('$y$ (mm)');
grid on;
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
% print(gcf,'figElastoPlasticCollapseEquivStr','-depsc');
% print(gcf,'figElastoPlasticCollapseEquivStr','-dpng');
% 
