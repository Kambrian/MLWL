figure;
plotgama_corr('H1','r.');
hold on;
plotgama_corr('H2','go');
plotgama_corr('H3','bs');
%%
figure;
plotgama_corr('H1.CC2GAMA','r.');
hold on;
plotgama_corr('H2.CC2GAMA','go');
plotgama_corr('H3.CC2GAMA','bs');
%%
figure;
plotgama_result('H1.CC2GAMA','n','r.');
hold on;
plotgama_result('H2.CC2GAMA','n','go');
plotgama_result('H3.CC2GAMA','n','bs');

%%
figure;
plotgama_corr('H14','.-')
figure;
plotgama_result('H14','n','.')
hold on;
plot(r,yiso{4}+1,'--');
plot(r,ynfw{4}+1,'-');
%% illustrate the necessity of averaging-correction
datadir='/mnt/charon/Lensing/output/';
name='H14';
file=[datadir,'shearprof_',name,'.txt'];
shear=importdata(file,'\t',1);
file=[datadir,'shearcov_',name,'.txt'];
cov=load(file);
file=[datadir,'randprof_',name,'.txt'];
rand=importdata(file,'\t',1);    
file=[datadir,'randcov_',name,'.txt'];
rcov=load(file);    
file=[datadir,'randcov_sw_',name,'.txt'];
rswcov=load(file);   
figure;
errorbar(shear.data(:,5),rand.data(:,4).*rand.data(:,end)./max(rand.data(:,end)),rand.data(:,8)./sqrt(rand.data(:,end)))
hold on; errorbar(shear.data(:,5),rand.data(:,4),rand.data(:,8)./sqrt(rand.data(:,end)),'r.')
set(gca,'xscale','log','yscale','log')