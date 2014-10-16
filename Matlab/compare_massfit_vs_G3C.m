load_G3Cv4
grp.MassProxyRaw=grp.MassProxy./(Ac+An./sqrt(grp.Mult)+Az./sqrt(grp.Zfof));
%%
grp.TotFluxRaw=grp.TotFluxProxy./(Bc+Bn./sqrt(grp.Mult)+Bz./sqrt(grp.Zfof));
%%
% newpar=[1.6,0,4.5];
% newpar=[4.15,-25.19,8.36];
newpar=[1500,-2500,42.6];
%Mnew=@(newpar) grp.MassProxyRaw.*(newpar(1)+newpar(2)./sqrt(grp.Mult)+newpar(3)./sqrt(grp.Zfof));
 Mnew=@(newpar) grp.TotFluxRaw.*(newpar(1)+newpar(2)./sqrt(grp.Mult)+newpar(3)./sqrt(grp.Zfof));
grp.NewMass=Mnew(newpar);
%%
flt=grp.Mult>2;
%%
figure;
plot(grp.LumMass,grp.LumMass,'r.','markersize',6);
hold on;
plot(grp.LumMass(flt),grp.NewMass(flt),'g.','markersize',6);

set(gca,'xscale','log','yscale','log')
xlabel('LumMass');ylabel('dynmass');legend('old','new')
% print('-depsc','../code/v5.0/python/comparenewmass_mult.eps')
%%
figure;
plot(grp.Mult,grp.LumMass,'r.','markersize',6);
hold on;
plot(grp.Mult,grp.NewMass,'g.','markersize',6);

% set(gca,'xscale','log','yscale','log')
xlabel('Mult');ylabel('dynmass');legend('old','new')
% print('-depsc','../code/v5.0/python/comparenewmass_mult.eps')
%%
figure;
plot(grp.Zfof,grp.LumMass,'r.','markersize',6);
hold on;
plot(grp.Zfof,grp.NewMass,'g.','markersize',6);
set(gca,'xscale','log','yscale','log')
xlabel('z');ylabel('dynmass');legend('old','new')
% print('-depsc','../code/v5.0/python/comparenewmass_z.eps')
%%
fold=@(n,z) (Ac+An./sqrt(n)+Az./sqrt(z));
fnew=@(p,n,z) (p(1)+p(2)./sqrt(n)+p(3)./sqrt(z));

z=0.1:0.05:0.5
for n=2:2:50
    plot(z,fold(n,z),'-')
    hold on;
    plot(z,fnew(newpar,n,z),'r:')
end
xlabel('z')
ylabel('Normalization factor')
%%
figure;
n=histc(grp.NewMass(flt)./grp.LumMass(flt),logspace(-1,1))
semilogx(logspace(-1,1),n)
