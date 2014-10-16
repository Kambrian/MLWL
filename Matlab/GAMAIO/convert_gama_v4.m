% cols=[2,3,4,6,7,8,9];%gamaid,ra,dec,zph,Q,zspec,class
% formats={[sprintf('\"'),'%d32',sprintf('\" ')],' %f ',' %f ','%f','%d','%f','%d'};
% file='/home/kam/Projects/Lensing/data/catgama_v4.csv';
% loadcsv_cols(file,cols,formats);

% save catgama_v4.mat CATA_INDEX   Q  SURVEY_CLASS  Z_SPEC   DEC  RA    Z_PHOTO

dir='/home/kam/Projects/Lensing/';
load catgama_v4 Z_SPEC Z_PHOTO Q SURVEY_CLASS

myfigure;
f=Q>2;
plot(Z_SPEC(f),Z_PHOTO(f),'.','markersize',3);
axis([0,1,0,1]);
hold on;
plot([0,1],[0,1],'r');
xlabel('$z_{GAMA}$');ylabel('SDSS $z_{ph}$');
print('-depsc',[dir,'output/zSpec_Phot_comp_gama4.eps']);

x=0:0.02:1;
myfigure;
[xm,ym,dyx]=linhist(Z_SPEC(f),x,'stairs','r');
hold on;
[xm,ym,dyx]=linhist(Z_PHOTO(f),x,'stairs');
xlabel('z');ylabel('Counts');
legend('Spec','Photo');
print('-depsc',[dir,'output/zSpec_Phot_distr_gama4.eps']);

%%
f=Q>2&Z_SPEC<1;
nbin=80;
myfigure;
[xx,yy,n,s]=densitygrid(Z_SPEC(f),Z_PHOTO(f),[nbin,nbin]);
contourf(xx,yy,log10((n+1)),10)
hold on;
plot([0,1],[0,1],'w');
axis([0,.8,0,.8]);
% imagesc(minmax(Z_SPEC(f)'),minmax(Z_PHOTO(f)'),log10(n+0.1));set(gca,'ydir','normal');
% colorbar;
% l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$z_{GAMA}$');ylabel('SDSS $z_{ph}$');
print('-depsc',[dir,'output/zSpec_Phot_distr2_gama4.eps']);