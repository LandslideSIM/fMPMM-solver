close all;clear;
set(0,'defaulttextinterpreter','latex') 
addpath('.\iterative_GIMP_solution\data','.\vectorized_GIMP_solution\data');   

for s = 1:6
vec=load(['GIMP_time_vectorized_' num2str(s) '.mat'])
loo=load(['GIMP_time_iterative_' num2str(s) '.mat'])

vRaw=vec.cycle_time;
vect=mean(vRaw);
lRaw=loo.cycle_time;
loop=mean(lRaw);
vmt = mean(vRaw(:,end-1)); vstd = std(vRaw(:,end-1));
lmt = mean(lRaw(:,end-1)); lstd = std(lRaw(:,end-1));


dv = vRaw(:,1:end-2)./repmat(mean(vRaw(:,1:end-2),1),length(vRaw),1);
dl = lRaw(:,1:end-2)./repmat(mean(lRaw(:,1:end-2),1),length(lRaw),1);


y = [loop(1:end-2);vect(1:end-2)];
yn = y ./ repmat([lmt;vmt],1,5)

fig4=figure
dim = [200 200 534 289];
set(fig4,'Units','pixels','Position',dim);
c = categorical({'iterative','vectorized'});
t = {'\texttt{SdS.m}','\texttt{p2Nsolve.m}','\texttt{mapN2p.m}','\texttt{DefUpdate.m}','\texttt{constitutive.m}'};
bar(c,y,'stacked');
strloop = ['$\bar{t}_c=',num2str(lmt,'%.2f'),'$ (s)'];
t1=text(0.80,1+lmt,strloop,'Interpreter','latex','Units','data','FontSize',14);
strvect = ['$\bar{t}_c=',num2str(vmt,'%.2f'),'$ (s)'];
t2=text(1.80,1+vmt,strvect,'Interpreter','latex','Units','data','FontSize',14);
ylabel('$\bar{t}_c$ (s)');
ylim([0 20])
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
h1=legend(t);
set(h1,'Interpreter','latex','Units','normalized','FontSize',12,'Position',[0.5299 0.5497 0.3613 0.3529],'NumColumns',1);
end

print(gcf,'figtimepercycle','-depsc');
print(gcf,'figtimepercycle','-dpng');



