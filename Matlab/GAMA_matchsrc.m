dir='/home/kam/Projects/Lensing/data/';
g{1}=loadGAMAcsv([dir,'groupsv2/cutcatG09.csv']);
g{2}=loadGAMAcsv([dir,'groupsv2/cutcatG12.csv']);
g{3}=loadGAMAcsv([dir,'groupsv2/cutcatG15.csv']);
gal=structcat([g{1};g{2};g{3}]);
% clear g;

% load_GAMACore;
load([dir,'GamaCoreDR1_v1.mat'],'GAMA_ID','SDSS_ID','RA_J2000','DEC_J2000','Z_HELIO','Z_QUALITY','Z_SOURCE','r_PETRO');
xgama=[RA_J2000,DEC_J2000];
zgama=[Z_HELIO,Z_QUALITY,Z_SOURCE];

%src
cd(dir);
% e=cell(3,1);
% e{1}=loadShearTBL(1,[1,2,4:8,16,53,57,58,59,60,54]);% ra,dec,run,rerun,camcol,field,id,r,z,z-sig,z+sig,z-2sig,z+2sig,z_template
% e{2}=loadShearTBL(2,[1,2,4:8,16,53,57,58,59,60,54]);
% e{3}=loadShearTBL(3,[1,2,4:8,16,53,57,58,59,60,54]);
% save source_id.mat e
load source_id.mat
%%
% %sdss
% sdss=cell(3,1);
% sdss{1}=fits_load_bintable('GAMAsrc1_kam.fit',[1:10,37:47]);
% sdss{2}=fits_load_bintable('GAMAsrc2_kam.fit',[1:10,37:47]);
% sdss{3}=fits_load_bintable('GAMAsrc3_kam.fit',[1:10,37:47]);
% save sdss_id.mat sdss
load sdss_id.mat
%%
% 1-percent or less galaxies have spectra
sky=3;
% spectral fraction
ngal=numel(sdss{sky}.photoz);
sum(~isnan(sdss{sky}.zsp))/ngal
% photoz success fraction, >0.99
sum(sdss{sky}.photoz>-100)/ngal
% zcc2 fraction, ~80 to 90 percent
sum(~isnan(sdss{sky}.photozcc2))/ngal
% zd1 fraction, ~80 to 90 percent
sum(~isnan(sdss{sky}.photozd1))/ngal

f=sdss{sky}.photoz>-100;
figure;plot(sdss{sky}.zsp(f),sdss{sky}.photoz(f),'.');
axis equal
%%
ll=cell(3,1);
for sky=1:3
    ll{sky}.ngrid=[40,20];
    [ll{sky}.grids,ll{sky}.xrange,ll{sky}.yrange,ll{sky}.step]=linklist(sdss{sky}.ra,sdss{sky}.dec,ll{sky}.ngrid);
end
%% match source to sdss check for IDs; 97.13% can be matched; 0.53% have multiple matches.
eps=10^-3.2;
nmatch=cell(3,1);
idmatch=cell(3,1);
dmin=cell(3,1);
ngal=zeros(3,1);
for sky=1:3
ngal(sky)=size(e{sky},1);
nmatch{sky}=zeros(ngal(sky),1);
idmatch{sky}=zeros(ngal(sky),1);
dmin{sky}=zeros(ngal(sky),1);
for i=1:ngal(sky)
    sources=select_grids(e{sky}(i,1:2),eps,[ll{sky}.xrange(1),ll{sky}.yrange(1)],ll{sky}.step,ll{sky}.ngrid,ll{sky}.grids);
    sources(abs(sdss{sky}.ra(sources)-e{sky}(i,1))>eps|abs(sdss{sky}.dec(sources)-e{sky}(i,2))>eps)=[];
    nmatch{sky}(i)=numel(sources);
    if nmatch{sky}(i)==1
       idmatch{sky}(i)=sources; 
       d=ccdist(e{sky}(i,1:2),[sdss{sky}.ra(sources),sdss{sky}.dec(sources)]);
       dmin{sky}(i)=acosd(d);
    elseif nmatch{sky}(i)>1
    d=ccdist(e{sky}(i,1:2),[sdss{sky}.ra(sources),sdss{sky}.dec(sources)]);
    [a,b]=min(d);
    idmatch{sky}(i)=sources(b);
    dmin{sky}(i)=acosd(a);
    end
end
end
% statistics
n=[nmatch{1};nmatch{2};nmatch{3}];
sum(n>0)/sum(ngal) %matchfrac
sum(n>1)/sum(ngal) % multi-match frac
d=[dmin{1};dmin{2};dmin{3}]; 
hist(log10(d(d>0))) % precision distribution

zmatch=cell(3,1);
zerrmatch=cell(3,1);
for sky=1:3
    zmatch{sky}=[sdss{sky}.photoz(idmatch{sky}(nmatch{sky}>0)),...
        sdss{sky}.photozcc2(idmatch{sky}(nmatch{sky}>0)),...
        sdss{sky}.photozd1(idmatch{sky}(nmatch{sky}>0))];
    zerrmatch{sky}=[sdss{sky}.photozerr(idmatch{sky}(nmatch{sky}>0)),...
        sdss{sky}.photozerrcc2(idmatch{sky}(nmatch{sky}>0)),...
        sdss{sky}.photozerrd1(idmatch{sky}(nmatch{sky}>0))];
end
% save src_sdss_match.mat nmatch idmatch dmin eps ngal zmatch zerrmatch
%% compare photoz's 
ze=[];zs=[];zsc=[];zsd=[];
for sky=1:3
    ze=[ze;e{sky}(nmatch{sky}>0,9)]; %zebra
    zs=[zs;sdss{sky}.photoz(idmatch{sky}(nmatch{sky}>0))]; %sdss template
    zsc=[zsc;sdss{sky}.photozcc2(idmatch{sky}(nmatch{sky}>0))]; %sdss cc2
    zsd=[zsd;sdss{sky}.photozd1(idmatch{sky}(nmatch{sky}>0))]; %sdss d1
end
myfigure;
[xx,yy,n,s]=densitygrid(ze(zs>-100),zs(zs>-100),[80,80]);
contourf(xx,yy,log10(n+1),20);
hold on;
plot([0,1],[0,1],'w');
xlabel('$z_{ZEBRA}$');ylabel('$z_{SDSS}$');
colorbar;
print('-depsc',[dir,'srccat-gama/photoz_zebra_sdss.eps']);
% ylim([0,1.5])

myfigure;
[xx,yy,n,s]=densitygrid(ze(zsc>-100),zsc(zsc>-100),[30,30]);
contourf(xx,yy,n);
hold on;
plot([0,1],[0,1],'w');
xlabel('$z_{ZEBRA}$');ylabel('$z_{SDSSCC2}$');
colorbar;
print('-depsc',[dir,'srccat-gama/photoz_zebra_sdssc.eps']);

myfigure;
[xx,yy,n,s]=densitygrid(ze(zsd>-100),zsd(zsd>-100),[30,30]);
contourf(xx,yy,n);
hold on;
plot([0,1],[0,1],'w');
xlabel('$z_{ZEBRA}$');ylabel('$z_{SDSSD1}$');
colorbar;
print('-depsc',[dir,'srccat-gama/photoz_zebra_sdssd.eps']);
%%  only half match for run and field, 70% match for camcol
ercf=cell(3,1); % run, camcol,field for e{}
srcf=cell(3,1); % run, camcol, field for sdss{}
for sky=1:3
ercf{sky}=e{sky}(logical(nmatch{sky}),[3,5,6]);
ids=idmatch{sky}(nmatch{sky}>0);
srcf{sky}=[sdss{sky}.run(ids),sdss{sky}.camcol(ids),sdss{sky}.field(ids)];
end
ercf=[ercf{1};ercf{2};ercf{3}];
srcf=[srcf{1};srcf{2};srcf{3}];

%% match gama gals to sdss: 99.8% can be matched; 90% have ANN photoz (cc2 and d1)
eps=10^-3.2;
ngal=numel(gal.CATA_INDEX);
nmatch=zeros(ngal,1);
idmatch=zeros(ngal,2);
dmin=nmatch;
for i=1:ngal
    if gal.RA(i)>120&&gal.RA(i)<150
        sky=1;
    elseif gal.RA(i)>150&&gal.RA(i)<200
            sky=2;
    elseif gal.RA(i)>200&&gal.RA(i)<240
                sky=3;
            else error('coordinate unexpected');
    end
    sources=select_grids([gal.RA(i),gal.DEC(i)],eps,[ll{sky}.xrange(1),ll{sky}.yrange(1)],ll{sky}.step,ll{sky}.ngrid,ll{sky}.grids);
    sources(abs(sdss{sky}.ra(sources)-gal.RA(i))>eps|abs(sdss{sky}.dec(sources)-gal.DEC(i))>eps)=[];
    nmatch(i)=numel(sources);
    if nmatch(i)==1
       idmatch(i,:)=[sources,sky]; 
       d=ccdist([gal.RA(i),gal.DEC(i)],[sdss{sky}.ra(sources),sdss{sky}.dec(sources)]);
       dmin(i)=acosd(d);
    elseif nmatch(i)>1
       d=ccdist([gal.RA(i),gal.DEC(i)],[sdss{sky}.ra(sources),sdss{sky}.dec(sources)]);
    [a,b]=min(d);
    idmatch(i,:)=[sources(b),sky];
    dmin(i)=acosd(a);
    end
end

save gama_sdss_match.mat nmatch idmatch dmin eps ngal
%% load match result
load gama_sdss_match.mat

sum(nmatch>1)/ngal
sum(nmatch>0)/ngal
figure;hist(log10(dmin(dmin>0)),20)
%% extract z and r
f=logical(nmatch);
j=0;
zsource=zeros(sum(f),6);
ztemp=zeros(sum(f),1)-1;
rsource=ztemp;
for i=1:ngal
    if nmatch(i)
        j=j+1;
        zsource(j,:)=[sdss{idmatch(i,2)}.photoz(idmatch(i,1)),sdss{idmatch(i,2)}.photoz(idmatch(i,1)),...
                                sdss{idmatch(i,2)}.photozcc2(idmatch(i,1)),sdss{idmatch(i,2)}.photozerrcc2(idmatch(i,1)),...
                                sdss{idmatch(i,2)}.photozd1(idmatch(i,1)),sdss{idmatch(i,2)}.photozerrd1(idmatch(i,1))];
        rsource(j)=sdss{idmatch(i,2)}.dered_r(idmatch(i,1));
        ztemp(j)=sdss{idmatch(i,2)}.pztype(idmatch(i,1),end);
    end
end
source_idmatch=cell(3,1);gama_zmatch=cell(3,1);
for i=1:3
source_idmatch{i}=idmatch(idmatch(:,2)==i,1);
gama_zmatch{i}=[gal.Z_SPEC(idmatch(:,2)==i),gal.SigErr(idmatch(:,2)==i)];
end
%stats
sum(zsource(:,3)>-100)/ngal
sum(zsource(:,5)>-100)/ngal
sum(zsource(:,1)>-100)/ngal
%% gama_v4 comparison: gama_photoz is not any of the sdss photoz's
gama4=load('catgama_v4','-mat','Z_SPEC','Z_PHOTO','Q','SURVEY_CLASS','CATA_INDEX');
id2ind=1:max(gama4.CATA_INDEX);
id2ind(gama4.CATA_INDEX)=1:numel(gama4.CATA_INDEX);

% load gama_sdss_match.mat

inds=id2ind(gal.CATA_INDEX(nmatch>0));
zsp4=gama4.Z_SPEC(inds);
zph4=gama4.Z_PHOTO(inds);
Q4=gama4.Q(inds);
CLASS4=gama4.SURVEY_CLASS(inds);

figure;
plot(gal.Z_SPEC(nmatch>0),zsp4,'.');
myfigure;
plot(zph4,zsource(:,1),'.');
myfigure;
plot(zph4,zsource(:,3),'.');
myfigure;
plot(zph4,zsource(:,5),'.');

ztmpl=[sdss{1}.photoz;sdss{2}.photoz;sdss{3}.photoz];
zcc2=[sdss{1}.photozcc2;sdss{2}.photozcc2;sdss{3}.photozcc2];
zd1=[sdss{1}.photozd1;sdss{2}.photozd1;sdss{3}.photozd1];
zsdss=zsource;  % matched part of sdss photozs
% save gama_redshifts.mat zsp4 zph4 zsdss ztmpl zcc2 zd1
%% fair comparison for the 4 photozs
v=log10(2):0.3:2.5;
f=zsdss(:,1)>-100&zsdss(:,3)>-100&zsdss(:,5)>-100;
z1=zsp4(f);
z2=zsdss(f,1);
myfigure;
ax1=axes('position',[0.2,0.5,0.3,0.3]);
[xx,yy,n,s]=densitygrid(z1,z2,[80,80],[0,0.5],[0,1]);
contourf(xx,yy,log10(n+1),v);
hold on;
plot([0,1],[0,1],'w');
% xlabel('$z_{GAMA}$');
ylabel('$z_{SDSStmpl}$');
axis([0,0.5,0,0.6]);

z2=zsdss(f,3);
ax2=axes('position',[0.5,0.5,0.3,0.3]);
[xx,yy,n,s]=densitygrid(z1,z2,[80,80],[0,0.5],[0,1]);
contourf(xx,yy,log10(n+1),v);
hold on;
plot([0,1],[0,1],'w');
% xlabel('$z_{GAMA}$');
ylabel('$z_{SDSScc2}$');
axis([0,0.5,0,0.6]);

z2=zsdss(f,5);
ax3=axes('position',[0.2,0.2,0.3,0.3]);
[xx,yy,n,s]=densitygrid(z1,z2,[80,80],[0,0.5],[0,1]);
contourf(xx,yy,log10(n+1),v);
hold on;
plot([0,1],[0,1],'w');
xlabel('$z_{GAMA}$');
ylabel('$z_{SDSSd1}$');
axis([0,0.5,0,0.6]);

z2=zph4(f);
ax4=axes('position',[0.5,0.2,0.3,0.3]);
[xx,yy,n,s]=densitygrid(z1,z2,[80,80],[0,0.5],[0,1]);
contourf(xx,yy,log10(n+1),v);
hold on;
plot([0,1],[0,1],'w');
xlabel('$z_{GAMA}$');
ylabel('$z_{GAMAph}$');
axis([0,0.5,0,0.6]);
set([ax1,ax2,ax3,ax4],'xtick',0:0.2:0.5,'xticklabel',0:0.2:0.4,'yticklabel',0:0.2:0.6);
set([ax1,ax2],'xticklabel','');
set([ax2,ax4],'yaxislocation','right');
print('-depsc',[dir,'srccat-gama/zgama_zphoto_conts.eps']);

%%
x=-1:0.01:1;
plttype='';
dz=[zsdss(f,1:2:5),zph4(f)]-repmat(zsp4(f),1,4);
figure;
linhist(dz(:,1),x,plttype,'r.-');hold on;
linhist(dz(:,2),x,plttype,'g.-');
linhist(dz(:,3),x,plttype,'bo-');
linhist(dz(:,4),x,plttype,'ks-');
figure;
h=cdfplot(-abs(dz(:,2)));
set(h,'color','r');
hold on;
h=cdfplot(-abs(dz(:,4)))
set(h,'color','b')
%% confidence regions
x=0:0.05:0.5;
alpha=0.683;
myfigure;
% [xm,ym,ylim]=skeleton(zsp4(f),zsdss(f,1),x,alpha);
% plot(xm,ym,'g-');hold on;
% plot(xm,ylim(:,1),'g--',xm,ylim(:,2),'g--');
[xm,ym,ylim]=skeleton(zsp4(f),zsdss(f,3),x,alpha);
h1=plot(xm,ym,'r-','linewidth',3);hold on;
plot(xm,ylim(:,1),'r--','linewidth',2);
plot(xm,ylim(:,2),'r--','linewidth',2);
% [xm,ym,ylim]=skeleton(zsp4(f),zsdss(f,5),x,alpha);
% plot(xm,ym,'b-');hold on;
% plot(xm,ylim(:,1),'b--',xm,ylim(:,2),'b--');
[xm,ym,ylim]=skeleton(zsp4(f),zph4(f),x,alpha);
h2=plot(xm,ym,'b-','linewidth',3);hold on;
plot(xm,ylim(:,1),'b--','linewidth',2);
plot(xm,ylim(:,2),'b--','linewidth',2);
plot([0,0.5],[0,0.5],'k-')    
axis([0,0.5,0,0.5]);

alpha=0.954;
[xm,ym,ylim]=skeleton(zsp4(f),zsdss(f,3),x,alpha);
plot(xm,ylim(:,1),'r:',xm,ylim(:,2),'r:');
[xm,ym,ylim]=skeleton(zsp4(f),zph4(f),x,alpha);
plot(xm,ylim(:,1),'b:',xm,ylim(:,2),'b:');
l=legend([h1,h2],'CC2','GAMA');set(l,'location','northwest');
xlabel('$z_{spec}$');ylabel('$z_{ph}$');
print('-depsc',[dir,'srccat-gama/PhotozCmp_Conf.eps']);
% 
% 
% figure;
% [xmed,ymed,ylim,xm,ym,ysig]=skeleton(zsp4(f),zsdss(f,3),x,alpha);
% errorbar(xm,ym,ysig,'ro-');
% hold on;
% [xmed,ymed,ylim,xm,ym,ysig]=skeleton(zsp4(f),zph4(f),x,alpha);
% errorbar(xm,ym,ysig,'bs-');
% plot([0,0.5],[0,0.5],'k:')    
%% joint distributions and dN/dz comparison
z1=gal.Z_SPEC(nmatch>0);
z2=zsource(:,1);
myfigure;
[xx,yy,n,s]=densitygrid(z1(z2>-100),z2(z2>-100),[80,80]);
contourf(xx,yy,log10(n+1),10);
hold on;
plot([0,1],[0,1],'w');
xlabel('$z_{GAMA}$');ylabel('$z_{SDSStmpl}$');
axis([min(xx(:)),0.5,min(yy(:)),0.8]);
% print('-depsc',[dir,'srccat-gama/zgama_zsdss_cont.eps']);
z2=zsource(:,3);
myfigure;
[xx,yy,n,s]=densitygrid(z1(z2>-100),z2(z2>-100),[80,80]);
contourf(xx,yy,log10(n+1),10);
hold on;
plot([0,1],[0,1],'w');
xlabel('$z_{GAMA}$');ylabel('$z_{SDSScc2}$');
axis([min(xx(:)),0.5,min(yy(:)),0.8]);
% print('-depsc',[dir,'srccat-gama/zgama_zsdsscc2_cont.eps']);
z2=zsource(:,5);
myfigure;
[xx,yy,n,s]=densitygrid(z1(z2>-100),z2(z2>-100),[80,80]);
contourf(xx,yy,log10(n+1),10);
hold on;
plot([0,1],[0,1],'w');
xlabel('$z_{GAMA}$');ylabel('$z_{SDSSzd1}$');
axis([min(xx(:)),0.5,min(yy(:)),0.8]);
print('-depsc',[dir,'srccat-gama/zgama_zsdsszd1_cont.eps']);

myfigure;
plot(gal.Z_SPEC(f),zsource(:,1),'marker','.','markersize',3,'linestyle','none');
axis([0,1,0,1]);
hold on;
plot([0,1],[0,1],'r-');
xlabel('$z_{GAMA}$');ylabel('$z_{SDSS}$');

x=0:0.02:1;
myfigure;
[xm,ym,dyx]=linhist(zsource(:,1),x);%,'stairsnorm','r-');
plot(xm,dyx/sum(ym),'r-');
hold on;
[xm,ym,dyx]=linhist(zsource(:,3),x);%,'stairsnorm','r-');
plot(xm,dyx/sum(ym),'g-');
[xm,ym,dyx]=linhist(zsource(:,5),x);%,'stairsnorm','r-');
plot(xm,dyx/sum(ym),'b-');
[xm,ym,dyx]=linhist(zph4,x);%,'stairsnorm','r-');
plot(xm,dyx/sum(ym),'c-');
plot(x+0.01,source_zdistr(x),'m:');
[xm,ym,dyx]=linhist(gal.Z_SPEC(f),x);%,'stairsnorm','g-');
plot(xm,dyx/sum(ym),'k-');
[xm,ym,dyx]=linhist([sdss{1}.photoz;sdss{2}.photoz;sdss{3}.photoz],x);%,'stairsnorm','r--');
plot(xm,dyx/sum(ym),'r--');
[xm,ym,dyx]=linhist([sdss{1}.photozcc2;sdss{2}.photozcc2;sdss{3}.photozcc2],x);%,'stairsnorm','r--');
plot(xm,dyx/sum(ym),'g--');
[xm,ym,dyx]=linhist([sdss{1}.photozd1;sdss{2}.photozd1;sdss{3}.photozd1],x);%,'stairsnorm','r--');
plot(xm,dyx/sum(ym),'b--');
[xm,ym,dyx]=linhist(gal.Z_SPEC,x);%,'stairsnorm','g--');
plot(xm,dyx/sum(ym),'k--');
l=legend('sdss tmpl','sdss cc2','sdss d1','gama photoz','M05','GAMA spec','All');
set(l,'box','off');
xlabel('z');
ylabel('dP/dz');
xlim([0,1]);
% print('-depsc',[dir,'srccat-gama/zgama_zsdss_distr.eps']);
%% magnitude dependence of z comparison
r=gal.PETROMAG_R_SDSS;
zspec=gal.Z_SPEC(f);
myfigure;
plot(gal.Z_SPEC(f),zsource(:,1),'.','markersize',3);hold on;
ff=(r<=19);ff=ff(f);
plot(zspec(ff),zsource(ff,1),'g.','markersize',3);
axis([0,0.6,0,1]);
hold on;
ff=(r<18);ff=ff(f);
plot(zspec(ff),zsource(ff,1),'r.','markersize',3);
plot([0,0.6],[0,0.6],'r-');
xlabel('$z_{GAMA}$');ylabel('$z_{ZEBRA}$');
% print('-depsc','GAMA_SRC_matchZ.eps');
%% error distribution: comparision to normal distribution;
z1=gal.Z_SPEC(nmatch>0);
z2=zsource(:,1);z2err=zsource(:,2);
z1=z1(z2>-100);z2err=z2err(z2>-100);z2=z2(z2>-100);
zvar=(z2-z1)./z2err;
z1=gal.Z_SPEC(nmatch>0);
z2=zsource(:,3);z2err=zsource(:,4);
z1=z1(z2>-100);z2err=z2err(z2>-100);z2=z2(z2>-100);
zvarcc2=(z2-z1)./z2err;
z1=gal.Z_SPEC(nmatch>0);
z2=zsource(:,5);z2err=zsource(:,6);
z1=z1(z2>-100);z2err=z2err(z2>-100);z2=z2(z2>-100);
zvard1=(z2-z1)./z2err;
x=-2:0.1:2;
myfigure;
[xm,ym,dyx]=linhist(zvar,x);
plot(xm,dyx/sum(ym),'r.');
hold on;
% [mu,sigma]=normfit(zvar)
% plot(x,normpdf(x,mu,sigma),'r-');
[xm,ym,dyx]=linhist(zvarcc2,x);
plot(xm,dyx/sum(ym),'go');
% [mu,sigma]=normfit(zvarcc2)
% plot(x,normpdf(x,mu,sigma),'g-');
[xm,ym,dyx]=linhist(zvard1,x);
plot(xm,dyx/sum(ym),'bx');
% [mu,sigma]=normfit(zvard1)
% plot(x,normpdf(x,mu,sigma),'b-');
plot(x,normpdf(x,0,1),'k-');
xlabel('$(z_{ph}-z_{sp})/\sigma_z$');
ylabel('Probability')
l=legend('SDSS Template','CC2','D1','Standard Normal'); set(l,'box','off','color','none','location','northwest');
print('-depsc',[dir,'srccat-gama/sdss_error_dist.eps']);
%% r magnitude comparison
myfigure;
plot(gal.PETROMAG_R_SDSS(f)-gal.A_r(f),rsource,'.');
% axis([10,25,10,25]);
hold on;
plot([12,20],[12,20],'r-');
xlabel('$r_{pet}$');ylabel('$r_{mod}$');
% print('-depsc','GAMA_SRC_matchR.eps');
