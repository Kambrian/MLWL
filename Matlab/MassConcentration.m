powerlaw_err=@(x,lnAerr,Berr,corrAB) sqrt(lnAerr^2+(log(x)*Berr).^2+2*log(x)*corrAB*lnAerr*Berr);
powerlawmodel_err=@(m,perr,c) powerlaw_err(m/1e12,perr(1)/log10(exp(1)),perr(2),c); %error in ln(y) for y=10^p(1)*(m/1e12).^p(2)
powerlawmodel=@(m,p) 10^p(1)*(m/1e12).^p(2);

p0=[1.006,-0.081];
x=logspace(12,15)';
y0=powerlawmodel(x,p0);
myfigure;
%% fix A*L^a at observational best fit (all redshift evolution fixed according to theory)
p=[1.36,-0.43];
perr=[0.68,0.52];
pcov=-0.98;
y=powerlawmodel(x,p);
yerr=powerlawmodel_err(x,perr,pcov);
yl=exp([log(y)+yerr,log(y)-yerr]);
htmp=area(x,[yl(:,1),yl(:,2)-yl(:,1)]);hold on;
set(htmp(2),'facecolor',[0.9,1,0.9],'edgecolor','w');
set(htmp(1),'facecolor','w','edgecolor','w');
h2=loglog(x,y,'g','displayname','M(L) freeze');hold on;
% plot(x,yl,'g:'); %error bounds
%% A*L^a + M-c
p=[1.67,-0.66];
perr=[0.33,0.25];
pcov=-0.91;
y=powerlawmodel(x,p);
yerr=powerlawmodel_err(x,perr,pcov);
yl=exp([log(y)-yerr,log(y)+yerr]);
set(gca,'layer','top');
h1=loglog(x,y,'r','displayname','M(L) free');hold on;
plot(x,yl,'r:'); %error bounds
%%
hobs=[];
x=logspace(log10(2e14),log10(4e15),20);
y=powerlawmodel(x/800,[log10(7.55),-0.45]);
[m,c]=convert_NFWmc(x,y,0,2,0);
h=loglog(m,c,'mx','displayname','SchmidtAllen07');hobs=[hobs,h];
x=logspace(log10(6e12),log10(2e15),20);
y=powerlawmodel(x/100,[log10(9.0),-0.172]);
[m,c]=convert_NFWmc(x,y,0,2,0);
h=loglog(m,c,'m+','displayname','Buote07');hobs=[hobs,h];
x=logspace(log10(5e13),log10(4e15),20);
y=powerlawmodel(x/13,[log10(14.5),-0.15]);
[m,c]=convert_NFWmc(x,y,0,2,0);
h=loglog(m,c,'m^','displayname','Comerford07');hobs=[hobs,h];
x=logspace(log10(5e12),log10(5e14),20);
y=powerlawmodel(x/13,[log10(4.1),-0.12]);
[m,c]=convert_NFWmc(x,y,1,2,0);
h=loglog(m,c,'bs','displayname','Johnston07');hobs=[hobs,h];
x=logspace(log10(3e13),log10(6e14),20);
y=powerlawmodel(x/100,[log10(4.6),-0.13]);
h=loglog(x,y,'bo','displayname','Mandelbaum08');hobs=[hobs,h];
x=logspace(log10(2e14),log10(1.5e15),20);
y=powerlawmodel(x/100,[log10(8.45),-0.41]);
[m,c]=convert_NFWmc(x,y,0,2,0);
h=loglog(m,c,'bd','displayname','Okabe09');hobs=[hobs,h];
%%
x=logspace(12,15)';
h0=loglog(x,y0,'k','displayname','Duffy08'); hold on;
l=legend([h0,h1,h2,hobs]);set(l,'interpreter','latex','location','southwest','fontsize',18,'box','off');
xlabel('$M[M_\odot/h]$');
ylabel('$c$');
xscale('log');
yscale('log');
ylim([1e-1,2e2]);
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-c.eps');