clear;
clc;
load G3Cv4up8/G3Cv4up8.mat
%% trim gals
names=fieldnames(gal);
f=gal.Rpetro<=19.4;
f=f&(~(gal.RA<142&gal.DEC<-0.7));%trim g09 for stellar mass completeness.
for j=1:numel(names)
    gal.(names{j})=gal.(names{j})(f);
end
%%
f=gal.CentralSampleIter>0&gal.SMsps>3e10&gal.SMsps<4e10;
f=f&gal.FlagSFR>0&gal.SSFR>10^-1.5;
data1.RA=gal.RA(f);
data1.DEC=gal.DEC(f);
data1.z=gal.Zspec(f);

data2.RA=gal.RA;
data2.DEC=gal.DEC;
data2.z=gal.Zspec;
%% random samples
[g1,g2,g3]=split_mockcat(data1);
rnd=GAMArand([2e4,2e4,2e4],-15,[numel(g1.RA),numel(g2.RA),numel(g3.RA)]);
rand1.RA=[rnd{1}(:,1);rnd{2}(:,1);rnd{3}(:,1)];
rand1.DEC=[rnd{1}(:,2);rnd{2}(:,2);rnd{3}(:,2)];
rand1.z=[g1.z(rnd{1}(:,3));g2.z(rnd{2}(:,3));g3.z(rnd{3}(:,3))];
%%
[g1,g2,g3]=split_mockcat(data2);
rnd=GAMArand([2e5,2e5,2e5],-25,[numel(g1.RA),numel(g2.RA),numel(g3.RA)]);
rand2.RA=[rnd{1}(:,1);rnd{2}(:,1);rnd{3}(:,1)];
rand2.DEC=[rnd{1}(:,2);rnd{2}(:,2);rnd{3}(:,2)];
rand2.z=[g1.z(rnd{1}(:,3));g2.z(rnd{2}(:,3));g3.z(rnd{3}(:,3))];
%%
rbin=logspace(-2,0,15);
d1d2=pair_counts(data1,data2,rbin);
% r=r./d1d2;
d1r2=pair_counts(data1,rand2,rbin);
r1d2=pair_counts(rand1,data2,rbin);
[r1r2,r]=pair_counts(rand1,rand2,rbin);
r=r./r1r2;
%%
n1=numel(data1.RA);n2=numel(data2.RA);
nr1=numel(rand1.RA);nr2=numel(rand2.RA);
ksiH=d1d2.*r1r2./d1r2./r1d2-1;
ksiLS=(d1d2/n1/n2-d1r2/n1/nr2-r1d2/nr1/n2)./(r1r2/nr1/nr2)+1;
figure;
plot(r,ksiH,'ro-');
hold on;
plot(r,ksiLS,'gx-');
xscale('log');
yscale('log');
ksiDP=d1d2/n2./(d1r2/nr2)-1;
ksiDP2=(d1d2/n1)./(r1d2/nr1)-1;
plot(r,ksiDP,'bs-',r,ksiDP2,'cd-');
%%
% save ksi.mat high_comov full_comov
%%
% load ksi.mat
% myfigure;
% plot(full_comov.r,full_comov.ksiLS,'ro-');
% hold on;
% plot(high_comov.r,high_comov.ksiLS,'gs-');
% xscale('log');
% yscale('log');
% xlabel('comoving R[Mpc/h]');
% ylabel('$\xi_{cg}$');
% legend('All Central','Active Central');
% title('$M_\star=(1\sim2)\times 10^{10}\mathrm{M}_\odot/h^2$');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/GAMASpecCorrelation.eps');