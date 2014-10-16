%%
% cd /work/Projects/Lensing/outputv4/data/
cd /mnt/charon/Lensing/output_v5.0/
fmts={'rd','gs','bo','kx'}
files={['WL_L1.DZ3.hdf5'],['WL_L2.DZ3.hdf5'],['WL_L3.DZ3.hdf5'],['WL_L3.DZ0.hdf5']};
%-----------
h=[];
myfigure;
for i=1:4
file=files{i}
m=h5read(file,'/predict1/Mmean')
z=h5read(file,'/predict1/z')
r=h5read(file,'/shear/seperation');
n=h5read(file,'/shear/numpair');
n=double(n);
nr=h5read(file,'/rand/numpair');
nerr=h5read(file,'/rand/numpair_err');
% n=h5read(file,'/shear/weight');
% nr=h5read(file,'/rand/weight');
% nerr=h5read(file,'/rand/weight_err');
nmock=h5read(file,'/rand/numrand');
nmock=double(nmock);
rv=comoving_200b_radius(m,0.3)/1000;
rat=n./(nr.*nmock/max(nmock));
rat_err=nerr./nr;
htmp=ploterr(r/rv,rat-1,[],rat_err*2,fmts{i},'logx');set(htmp(1),'markersize',10);
h=[h;htmp(1)];
hold on;
end
l=legend(h,'L1','L2','L3','L3Tight'); set(l,'interpreter','latex');
plot([1e-2,30],[0,0],'k--')
xlabel('$R/R_{v}$');ylabel('$n/n_{rand}\frac{\ }{\ }1$');
set(gca,'xscale','log');
xlim([1e-2,20])
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/Contamination.eps')