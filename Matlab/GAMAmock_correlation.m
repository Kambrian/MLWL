clear
clc
load G3Cv4up8/mockcat_1.mat
%% trim gals
    raminmax=[129,141;
        174,186;
        211.5,223.5];
%     decminmax=[-1,3;-2,2;-2,2];
    decminmax=[-0.7,3;-2,2;-2,2];
names=fieldnames(galmock);
f=zeros(size(galmock.RA));
for sky=1:3
    f=f|(galmock.RA>raminmax(sky,1)&galmock.RA<=raminmax(sky,2)&galmock.DEC>decminmax(sky,1)&galmock.DEC<=decminmax(sky,2));
end
f=f&galmock.r_mag<=19.4;
for j=1:numel(names)
    galmock.(names{j})=galmock.(names{j})(f);
end
%%
f=galmock.CentralSampleIter>0&galmock.Mstar>3e10&galmock.Mstar<4e10;
% f=f&galmock.SSFR>10^-1.5;
data1.RA=galmock.RA(f);
data1.DEC=galmock.DEC(f);
data1.z=galmock.Z(f);

data2.RA=galmock.RA;
data2.DEC=galmock.DEC;
data2.z=galmock.Z;
%% random samples
[g1,g2,g3]=split_mockcat(data1);
rnd=GAMArand([2e4,2e4,2e4],-10,[numel(g1.RA),numel(g2.RA),numel(g3.RA)]);
rand1.RA=[rnd{1}(:,1);rnd{2}(:,1);rnd{3}(:,1)];
rand1.DEC=[rnd{1}(:,2);rnd{2}(:,2);rnd{3}(:,2)];
rand1.z=[g1.z(rnd{1}(:,3));g2.z(rnd{2}(:,3));g3.z(rnd{3}(:,3))];
%%
[g1,g2,g3]=split_mockcat(data2);
rnd=GAMArand([2e5,2e5,2e5],-20,[numel(g1.RA),numel(g2.RA),numel(g3.RA)]);
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
% mock_full_comov.r=r;
% mock_full_comov.ksiH=ksiH;
% mock_full_comov.ksiLS=ksiLS;
% mock_full_comov.ksiDP=ksiDP;
% save ksi.mat mock_full_comov -append
%%
% load ksi2.mat
myfigure;
plot(full_comov.r,full_comov.ksiLS,'ro-','markerfacecolor','r');
hold on;
plot(high_comov.r,high_comov.ksiLS,'gs-','markerfacecolor','g');
plot(mock_full_comov.r,mock_full_comov.ksiLS,'ro--');
hold on;
plot(mock_high_comov.r,mock_high_comov.ksiLS,'gd--');
xscale('log');
yscale('log');
xlabel('comoving R[Mpc/h]');
ylabel('$\xi_{cg}$');
l=legend('All Central','Active Central','All Central(Mock)','Active Central(Mock)');set(l,'interpreter','latex');
title('$M_\star=(3\sim4)\times 10^{10}\mathrm{M}_\odot/h^2,\Delta z<0.005$');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/GAMASpecCorrelation2.eps');