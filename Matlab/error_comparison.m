datadir='/mnt/charon/Lensing/output/';
name='H3';
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
file=[datadir,'bootprof_',name,'.txt'];
boot=importdata(file,'\t',1);    
file=[datadir,'bootcov_',name,'.txt'];
bootcov=load(file);    
file=[datadir,'bootcov_sw_',name,'.txt'];
bootswcov=load(file);   
file=[datadir,'sbootprof_',name,'.txt'];
sboot=importdata(file,'\t',1);    
file=[datadir,'sbootcov_',name,'.txt'];
sbootcov=load(file);    
file=[datadir,'sbootcov_sw_',name,'.txt'];
sbootswcov=load(file);   
file=[datadir,'sphotoprof_',name,'.txt'];
sphoto=importdata(file,'\t',1);    
file=[datadir,'sphotocov_',name,'.txt'];
sphotocov=load(file);    
file=[datadir,'sphotocov_sw_',name,'.txt'];
sphotoswcov=load(file);   
file=[datadir,'sshuffleprof_',name,'.txt'];
sshuffle=importdata(file,'\t',1);    
file=[datadir,'sshufflecov_',name,'.txt'];
sshufflecov=load(file);    
file=[datadir,'sshufflecov_sw_',name,'.txt'];
sshuffleswcov=load(file);   
file=[datadir,'srotateprof_',name,'.txt'];
srotate=importdata(file,'\t',1);    
file=[datadir,'srotatecov_',name,'.txt'];
srotatecov=load(file);    
file=[datadir,'srotatecov_sw_',name,'.txt'];
srotateswcov=load(file);   
file=[datadir,'bootchunkprof_',name,'.txt'];
bootchunk=importdata(file,'\t',1);    
file=[datadir,'bootchunkcov_',name,'.txt'];
bootchunkcov=load(file);    
file=[datadir,'bootchunkcov_sw_',name,'.txt'];
bootchunkswcov=load(file);   
file=[datadir,'sbootchunkprof_',name,'.txt'];
sbootchunk=importdata(file,'\t',1);    
file=[datadir,'sbootchunkcov_',name,'.txt'];
sbootchunkcov=load(file);    
file=[datadir,'sbootchunkcov_sw_',name,'.txt'];
sbootchunkswcov=load(file);   
%%
figure;
plot(shear.data(:,5),shear.data(:,4),'r.');
hold on;
plot(shear.data(:,5),sqrt(diag(cov)),'go');
plot(shear.data(:,5),sqrt(diag(rcov)),'bx');
plot(shear.data(:,5),sqrt(diag(srotatecov)),'ks');
set(gca,'xscale','log');
set(gca,'yscale','log');
%%
myfigure;
plot(shear.data(:,5),sqrt(diag(cov))./shear.data(:,4),'r.-');
hold on;
plot(shear.data(:,5),sqrt(diag(rcov))./shear.data(:,4),'ko-');
plot(shear.data(:,5),sqrt(diag(bootcov))./shear.data(:,4),'bo-');
plot(shear.data(:,5),sqrt(diag(sbootcov))./shear.data(:,4),'bs-');
plot(shear.data(:,5),sqrt(diag(srotatecov))./shear.data(:,4),'c*-');
plot(shear.data(:,5),sqrt(diag(sshufflecov))./shear.data(:,4),'cx-');
plot(shear.data(:,5),sqrt(diag(sphotocov))./shear.data(:,4),'m>-');
plot(shear.data(:,5),sqrt(diag(bootchunkcov))./shear.data(:,4),'gp-');
plot(shear.data(:,5),sqrt(diag(sbootchunkcov))./shear.data(:,4),'gd-');
l=legend('Analytical','RandLens','LensBoot','SrcBoot','ShearRotate','ShearShuffle','PhotoZ-MC','LensChunkBoot','SrcChunkBoot');
set(l,'location','northwest','box','off','fontsize',15);
xlabel('r/(Mpc/h)');ylabel('$\sigma/\sigma_{simple}$');
set(gca,'xscale','log');
% set(gca,'xlim',[0.01,300],'xtick',[0.01,0.1,1,10,100]);
print('-depsc',['error_comp_',name,'.eps']);

%%
myfigure;
plot(shear.data(:,5),shear.data(:,4)./sqrt(diag(cov)),'r.-');
hold on;
plot(shear.data(:,5),sqrt(diag(rcov))./sqrt(diag(cov)),'ko-');
plot(shear.data(:,5),sqrt(diag(bootcov))./sqrt(diag(cov)),'bo-');
plot(shear.data(:,5),sqrt(diag(sbootcov))./sqrt(diag(cov)),'bs-');
plot(shear.data(:,5),sqrt(diag(srotatecov))./sqrt(diag(cov)),'c*-');
plot(shear.data(:,5),sqrt(diag(sshufflecov))./sqrt(diag(cov)),'cx-');
plot(shear.data(:,5),sqrt(diag(sphotocov))./sqrt(diag(cov)),'m>-');
plot(shear.data(:,5),sqrt(diag(bootchunkcov))./sqrt(diag(cov)),'gp-');
plot(shear.data(:,5),sqrt(diag(sbootchunkcov))./sqrt(diag(cov)),'gd-');
l=legend('Simple','RandLens','LensBoot','SrcBoot','ShearRotate','ShearShuffle','PhotoZ-MC','LensChunkBoot','SrcChunkBoot');
% l=legend('Simple','RandLens','LensBoot','SrcBoot','ShearRotate','ShearShuffle','LensChunkBoot','SrcChunkBoot');

set(l,'location','southwest','box','off','fontsize',15);
xlabel('r/(Mpc/h)');ylabel('$\sigma/\sigma_{analytical}$');
set(gca,'xscale','log');
set(gca,'xlim',[0.01,300],'xtick',[0.01,0.1,1,10,100]);
print('-depsc',['error_comp2_',name,'.eps']);
%%
R=cov;
C=cov;
for i=1:size(C,1)
    for j=1:size(C,2)
        R(i,j)=C(i,j)/sqrt(C(i,i))/sqrt(C(j,j));
    end
end
figure; imagesc(R);colorbar();%colormap('gray');
print('-depsc',['cov_',name,'.eps']);

R=cov;
C=srotatecov;
for i=1:size(C,1)
    for j=1:size(C,2)
        R(i,j)=C(i,j)/sqrt(C(i,i))/sqrt(C(j,j));
    end
end
figure; imagesc(R);colorbar();%colormap('gray');
print('-depsc',['cov_srotate_',name,'.eps']);
figure;imagesc(cov./C);colorbar();

%%
figure;
loglog(shear.data(:,5),diag(cov).*shear.data(:,7),'.')