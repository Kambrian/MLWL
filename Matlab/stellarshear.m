myfigure;
l=importdata('shearprof_GALlow.txt','\t',1)
h=importdata('shearprof_GALhigh.txt','\t',1)
errorbar(l.data(:,4),l.data(:,1),l.data(:,3),'rs');
hold on;
errorbar(h.data(:,4),h.data(:,1),h.data(:,3),'ro');
sl=importdata('shearprof_satlow.txt','\t',1)
sh=importdata('shearprof_sathigh.txt','\t',1)
errorbar(sl.data(:,4),sl.data(:,1),sl.data(:,3),'gs');
hold on;
errorbar(sh.data(:,4),sh.data(:,1),sh.data(:,3),'go');
bl=importdata('shearprof_BCGlow.txt','\t',1)
bh=importdata('shearprof_BCGhigh.txt','\t',1)
errorbar(bl.data(:,4),bl.data(:,1),bl.data(:,3),'ks-');
hold on;
errorbar(bh.data(:,4),bh.data(:,1),bh.data(:,3),'ko-');
set(gca,'xscale','log','yscale','log');
legend('M_*<1.05e11','M_*>1.05e11');
xlabel('comoving r/(Mpc/h)');ylabel('$\Delta\Sigma$');
print('-depsc','wlSig_GAL_more.eps');
%%
myfigure;
f=importdata('shearprof_field.txt','\t',1)
errorbar(f.data(:,4),f.data(:,1),f.data(:,3),'ro');
s=importdata('shearprof_sat.txt','\t',1)
hold on;
errorbar(s.data(:,4),s.data(:,1),s.data(:,3),'go');
b=importdata('shearprof_BCG.txt','\t',1)
hold on;
errorbar(b.data(:,4),b.data(:,1),b.data(:,3),'ko-');
set(gca,'xscale','log','yscale','log');
legend('field','sat','bcg');
xlabel('comoving r/(Mpc/h)');ylabel('$\Delta\Sigma$');
print('-depsc','wlSig_GAL_env.eps');
%%

lr=importdata('randprof_GALlow.txt','\t',1)
hr=importdata('randprof_GALhigh.txt','\t',1)

myfigure;
err_l=l.data(:,6)./lr.data(:,3).*sqrt(1./l.data(:,6)+(lr.data(:,6)./sqrt(lr.data(:,7))./lr.data(:,3)).^2);
err_h=h.data(:,6)./hr.data(:,3).*sqrt(1./h.data(:,6)+(hr.data(:,6)./sqrt(hr.data(:,7))./hr.data(:,3)).^2);
errorbar(l.data(:,4),l.data(:,6)./lr.data(:,3),err_l,'rs');
hold on;
errorbar(h.data(:,4),h.data(:,6)./hr.data(:,3),err_h,'go');
set(gca,'xscale','log');
plot(get(gca,'xlim'),[1,1]);
legend('M_*<1.05e11','M_*>1.05e11');
xlabel('comoving r/(Mpc/h)');ylabel('$n/n_{rand}$');

% print('-depsc','nprof_rat_GAL.eps');

myfigure;
semilogx(l.data(:,4),lr.data(:,4)./l.data(:,3),'.');
hold on;
semilogx(h.data(:,4),hr.data(:,4)./h.data(:,3),'r.');
