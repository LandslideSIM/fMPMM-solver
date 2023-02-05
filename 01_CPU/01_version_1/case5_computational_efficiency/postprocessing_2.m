close all;clear;
set(0,'defaulttextinterpreter','latex') 
addpath('.\iterative_GIMP_solution\data','.\vectorized_GIMP_solution\data');   

for s = 1:6

vec=load(['GIMP_time_vectorized_' num2str(s) '.mat']);
loo=load(['GIMP_time_iterative_' num2str(s) '.mat']);

vectorized(s)=mean(vec.runt(:,1));
iterative(s)=mean(loo.runt(:,1));

vec=load(['GIMP_vectorized_data_' num2str(s) '.mat']);
loo=load(['GIMP_iterative_data_' num2str(s) '.mat']);

nelx(s)=vec.meD.nEx-3;
nely(s)=vec.meD.nEy-3;
np(s) = vec.mpD.n    ;

data_vectorized(s) = min(vec.meD.h);
data_iterative(s) = min(loo.meD.h);
end

fig1=figure(1)
set(fig1,'Units','pixels','Position',[200 200 541 277]);
hold on
ax1=plot(1./data_vectorized,vectorized,'bo--','LineWidth',2);
ax2=plot(1./data_iterative,iterative,'rs--','LineWidth',2);
hold off
set(gca, 'XScale', 'log', 'YScale', 'log')
box on
grid on
%xlim([100 150000])
ylim([4e-2 1e3])
tit = {'Vectorized','Iterative'};
h1=legend([ax1 ax2],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.6415 0.6440 0.2420 0.2383],'NumColumns',1);
title(h1,'$n_{pe}=3^2$');
yticks([1e0 1e2])
xlabel('$1/h$ (m)');
ylabel('it/s');
grid on;
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
print(gcf,'figItpsnp','-depsc');
print(gcf,'figItpsnp','-dpng');


