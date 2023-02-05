clear,close
set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;
nstr = 4;

for run = 1:7
D_vect{run} = load(['.\data_save\data_GIMPM\GIMP_time_vectorized_' num2str(run) '.mat']);
GIMP = load(['.\data_save\data_GIMPM\GIMP_data_' num2str(run) '.mat']);
np(run) = GIMP.mpD.n;
nNe     = GIMP.meD.nNe;
nN(run) = GIMP.meD.nN;
nelx(run)= GIMP.meD.nEx-4;
nely(run)= GIMP.meD.nEy-4;
tsolve_vect(run) = D_vect{run}.tsolve
cycle_vect{run} = D_vect{run}.cycle_time

tm(run,:,1)   = mean(cycle_vect{run},1)
FLOP(run,1,1) = (np(run)*nNe  )*11;
FLOP(run,2,1) = (np(run)*nNe  )*1 +(np(run)*nNe*2     )*2+(nN(run)*2    )*5+(nstr*nNe*2*np(run))*1;
FLOP(run,3,1) = (np(run)*nNe*2)*2 +(np(run)           )*2+(np(run)*nNe*2)*1+(nN(run)*2         )*1;
FLOP(run,4,1) = (np(run)*nNe  )*6 +(nstr*nNe*2*np(run))*1+(np(run)      )*11;
FLOP(run,5,1) = (np(run)      )*16+(nstr*np(run)      )*6;

if(run<7)
    D_iter{run} = load(['.\data_save\data_GIMPM_iterative\GIMP_time_vectorized_' num2str(run) '.mat']);
    tsolve_iter(run) = D_iter{run}.tsolve
    cycle_iter{run} = D_iter{run}.cycle_time
    
    tm(run,:,2)   = mean(cycle_iter{run},1)
 
    FLOP(run,1,2) = (np(run)*nNe  )*11;
    FLOP(run,2,2) = (np(run)*nNe  )*1 +(np(run)*nNe*2     )*2+(nN(run)*2    )*5+(nstr*nNe*2*np(run))*1;
    FLOP(run,3,2) = (np(run)*nNe*2)*2 +(np(run)           )*2+(np(run)*nNe*2)*1+(nN(run)*2         )*1;
    FLOP(run,4,2) = (np(run)*nNe  )*6 +(nstr*nNe*2*np(run))*1+(np(run)      )*11;
    FLOP(run,5,2) = (np(run)      )*16+(nstr*np(run)      )*6;
end
end

FLOPS = FLOP(:,1:5,:)./tm(:,1:5,:)
fig3=figure(3);
clf;
set(fig3,'Units','pixels','Position',[100 100 541 277]);
hold on;
ax1=loglog(np,sum(FLOPS(:,1,1),2)/1e6,'s--');
ax2=loglog(np,sum(FLOPS(:,2,1),2)/1e6,'d--');
ax3=loglog(np,sum(FLOPS(:,3,1),2)/1e6,'x--');
ax4=loglog(np,sum(FLOPS(:,4,1),2)/1e6,'*--');
ax5=loglog(np,sum(FLOPS(:,5,1),2)/1e6,'+--');
ax6=loglog(np,sum(FLOPS(:,1:5,1),2)/1e6,'ro-');
hold off;
box on;
xlim([7.5890 7.5890e+04])
ylim([1.7944 1.7944e+03])
yticks([1e1 1e2 1e3])
set(gca, 'YScale', 'log','XScale','log')
tit = {'\texttt{SdS.m}','\texttt{p2Nsolve.m}','\texttt{mapN2p.m}','\texttt{DefUpdate.m}','\texttt{constitutive.m}','solver'};
h1=legend([ax1 ax2 ax3 ax4 ax5 ax6],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.3340 0.2322 0.5615 0.2265],'NumColumns',2);
xlabel('$n_p$ (-)')
ylabel('Mflops')
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
print(gcf,'figElastoPlasticCollapseMflops','-depsc');

fig1=figure(1);
clf;
set(fig1,'Units','pixels','Position',[100 100 541 277]);
hold on;
ax1=plot(np         ,1./tm(:  ,6,1),'o-');
ax2=plot(np(1:end-1),1./tm(1:6,6,2),'s-');
% ax3=plot(np         ,ones(size(np))*1   ,'r:');
% ax3=plot(np         ,ones(size(np))*60  ,'r:');
% ax3=plot(np         ,ones(size(np))*3600,'r:');
hold off;
box on;
tit = {'Vectorised solver','Iterative solver'};
h1=legend([ax1 ax2],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.1495 0.2525 0.3802 0.1553],'NumColumns',1);
set(gca, 'YScale', 'log','XScale','log')
yticks([1 10 100 1000])
xlabel('$n_p$ (-)')
ylabel('it/s')
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
ylim([0.30 1500])
xlim([8 60000])
print(gcf,'figElastoPlasticCollapseItpersec','-depsc');


(1./tm(1:6,6,1))./(1./tm(1:6,6,2))











for i=1:7
Y_vect(i,:) = mean(cycle_vect{i},1)
if(i<7)
Y_iter(i,:) = mean(cycle_iter{i},1)
end
end
Y_vect = Y_vect(:,1:end-2);
Y_iter = Y_iter(:,1:end-2);
fig2=figure(2);
clf;
set(fig2,'Units','pixels','Position',[100 100 541 277]);
hold on
ax4=plot(np(1:end-1),1./Y_iter(:,1),'rs--')
ax5=plot(np(1:end-1),1./Y_iter(:,2),'bs--')
ax6=plot(np(1:end-1),1./Y_iter(:,3),'gs--')
ax1=plot(np         ,1./Y_vect(:,1),'ro-' )
ax2=plot(np         ,1./Y_vect(:,2),'bo-' )
ax3=plot(np         ,1./Y_vect(:,3),'go-' )
set(gca, 'YScale', 'log','XScale','log')
hold off
box on
axis tight
tit = {'\texttt{SdS}','\texttt{p2Nsolve}','\texttt{mapN2p}'};
h1=legend([ax1 ax2 ax3],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.6767 0.5954 0.2185 0.3097],'NumColumns',1);
title(h1,'Vectorised')
ah1 = axes('position',get(gca,'position'),'visible','off')
h2=legend(ah1,[ax4 ax5 ax6],tit);
set(h2,'Interpreter','latex','FontSize',12,'Position',[0.1405 0.1314 0.2185 0.3097],'NumColumns',1);
title(h2,'Iterative')
box on







