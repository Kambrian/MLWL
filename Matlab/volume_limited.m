declim=[-1,3;-2,2;-2,2];
Area=12*pi/180*sum(sind(declim(:,2))-sind(declim(:,1)));
zbin=[0:0.03:0.5]';
zmid=(zbin(1:end-1)+zbin(2:end))/2;
rbin=comoving_dist(0.3,0.7,zbin);
Vbin=Area/3*diff(rbin.^3);
%%
xbrd=[%5e9,1e10;
    1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 
    %3.7e11, 6e11
    ];
% c=['rgbcmykrgbcmyk'];
c=spring(size(xbrd,1)+1);
myfigure;
for i=1:size(xbrd,1)
f=gal.SMsps>xbrd(i,1)&gal.SMsps<xbrd(i,2)&gal.CentralSampleIter>0;
N=histc(gal.Zspec(f),zbin);
n=N(1:end-1)./Vbin; %comoving density, (Mpc/h)^-3
en=sqrt(N(1:end-1))./Vbin;
Ndens=[zmid,n,en];
% save GAMAI-Ndens-z.txt Ndens -ascii
plot(zmid,n,'color',c(i,:), 'displayname', ['$',printexp10(mean(xbrd(i,:))),'$']);hold on;
end
l=legend('show');set(l,'interpreter','latex');
yscale('log');
y=ylim();
plot([0.2,0.2],y,'k--');
xlabel('z');
ylabel('comoving density [(Mpc/$h)^{-3}$]'); 

print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/dNdz.eps');