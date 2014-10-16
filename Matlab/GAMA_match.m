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
%all gama galaxies are within GamaCore
nmax=max(max(GAMA_ID),max(gal.CATA_INDEX));
fl=zeros(nmax,1);
fr=fl;
fl(GAMA_ID)=1;
fr(gal.CATA_INDEX)=1;  
f=fl-fr;
find(f<0)
sum(f)
id2gama=zeros(nmax,1);
id2gama(GAMA_ID)=1:numel(GAMA_ID);
% zbad=Z_QUALITY(id2gama(logical(f)));

%%
ll=cell(3,1);
for sky=1:3
    ll{sky}.ngrid=[40,20];
    [ll{sky}.grids,ll{sky}.xrange,ll{sky}.yrange,ll{sky}.step]=linklist(e{sky}(:,1),e{sky}(:,2),ll{sky}.ngrid);
end
%% match GAMA_Core to source to check for IDs
eps=10^-3.2;
ngal=numel(GAMA_ID);
nmatch=zeros(ngal,1);
idmatch=zeros(ngal,2);
dmin=nmatch;
for i=1:ngal
    if RA_J2000(i)>120&&RA_J2000(i)<150
        sky=1;
    elseif RA_J2000(i)>150&&RA_J2000(i)<200
            sky=2;
    elseif RA_J2000(i)>200&&RA_J2000(i)<240
                sky=3;
            else error('coordinate unexpected');
    end
    sources=select_grids([RA_J2000(i),DEC_J2000(i)],eps,[ll{sky}.xrange(1),ll{sky}.yrange(1)],ll{sky}.step,ll{sky}.ngrid,ll{sky}.grids);
    sources(abs(e{sky}(sources,1)-RA_J2000(i))>eps|abs(e{sky}(sources,2)-DEC_J2000(i))>eps)=[];
    nmatch(i)=numel(sources);
    if nmatch(i)==1
       idmatch(i,:)=[sources,sky]; 
    elseif nmatch(i)>1
    d=ccdist([RA_J2000(i),DEC_J2000(i)],e{sky}(sources,[1,2]));
    [a,b]=min(d);
    idmatch(i,:)=[sources(b),sky];
    dmin(i)=acosd(a);
    end
end
sid=SDSS_ID(logical(nmatch));
j=0;
run=zeros(size(sid));field=run;camcol=run;
for i=1:ngal
    if nmatch(i)
        j=j+1;
        run(j)=e{idmatch(i,2)}(idmatch(i,1),3);
        camcol(j)=e{idmatch(i,2)}(idmatch(i,1),5);
        field(j)=e{idmatch(i,2)}(idmatch(i,1),6);
    end
end
runmask=hex2dec('0000FFFF00000000');
cammask=hex2dec('00000000E0000000');
fieldmask=hex2dec('000000000FFF0000');
grun=bitshift(bitand(sid,runmask),-32);
gcamcol=bitshift(bitand(sid,cammask),-29);
gfield=bitshift(bitand(sid,fieldmask),-16);
grun=double(grun);gcamcol=double(gcamcol);gfield=double(gfield);
%ObjID bit en-coding:
%http://cas.sdss.org/dr5/en/help/docs/algorithm.asp?search=objid&submit1=Search
%
%0 	1 	0x8000000000000000 	empty 	unassigned
%1-4 	4 	0x7800000000000000 	skyVersion 	resolved sky version (0=TARGET, 1=BEST, 2-15=RUNS)
%5-15 	11 	0x07FF000000000000 	rerun 	number of pipeline rerun
%16-31 	16 	0x0000FFFF00000000 	run 	run number
%32-34 	3 	0x00000000E0000000 	camcol 	camera column (1-6)
%35 	1 	0x0000000010000000 	firstField 	is this the first field in segment?
%36-47 	12 	0x000000000FFF0000 	field 	field number within run
%48-63 	16 	0x000000000000FFFF 	object 	object number within field

%% match gama gals to source 
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
    sources(abs(e{sky}(sources,1)-gal.RA(i))>eps|abs(e{sky}(sources,2)-gal.DEC(i))>eps)=[];
    nmatch(i)=numel(sources);
    if nmatch(i)==1
       idmatch(i,:)=[sources,sky]; 
       d=ccdist([gal.RA(i),gal.DEC(i)],e{sky}(sources,[1,2]));
       dmin(i)=acosd(d);
    elseif nmatch(i)>1
    d=ccdist([gal.RA(i),gal.DEC(i)],e{sky}(sources,[1,2]));
    [a,b]=min(d);
    idmatch(i,:)=[sources(b),sky];
    dmin(i)=acosd(a);
    end
end
sum(nmatch>1)/ngal
sum(nmatch>0)/ngal
max(dmin)
figure;hist(log10(dmin(dmin>0)),20)
f=logical(nmatch);
j=0;
zsource=zeros(sum(f),5);
ztemp=zeros(sum(f),1)-1;
rsource=ztemp;
for i=1:ngal
    if nmatch(i)
        j=j+1;
        zsource(j,:)=e{idmatch(i,2)}(idmatch(i,1),9:13);
        rsource(j)=e{idmatch(i,2)}(idmatch(i,1),8);
        ztemp(j)=e{idmatch(i,2)}(idmatch(i,1),end);
    end
end
source_idmatch=cell(3,1);gama_zmatch=cell(3,1);
for i=1:3
source_idmatch{i}=idmatch(idmatch(:,2)==i,1);
gama_zmatch{i}=[gal.Z_SPEC(idmatch(:,2)==i),gal.SigErr(idmatch(:,2)==i)];
end
% save GAMA_SRC_matchZ.mat source_idmatch gama_zmatch
%% redshift distribution and comparison
myfigure;
% set(gcf,'paperunits','centimeters','paperposition',[0.6,5,20,17]);
ax1=axes('position',[0.1607    0.1119    0.7443    0.7631]);
line(gal.Z_SPEC(f),zsource(:,1),'marker','.','markersize',3,'linestyle','none');
axis([0,0.6,0,1]);
hold on;
plot([0,0.6],[0,0.6],'r-');
xlabel('$z_{GAMA}$');ylabel('$z_{ZEBRA}$');
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'Color','none',...
           'XColor','r');
[xm,ym,dyx]=linhist(zsource(:,1),40);
h1=line(dyx/sum(ym),xm,'color','r','linestyle','-');
hold on;
[xm,ym,dyx]=linhist([e{1}(:,9);e{2}(:,9);e{3}(:,9)],40);
h2=line(dyx/sum(ym),xm,'color','r','linestyle','--');
ylim([0,1]);
xlabel('$dP/dz_{ZEBRA}$');
ax3 = axes('Position',get(ax1,'Position'),...
           'YAxisLocation','right',...
           'Color','none',...
           'YColor','g');       
[xm,ym,dyx]=linhist(gal.Z_SPEC(f),30);  
line(xm,dyx/sum(ym),'color','g','linestyle','-');
hold on;
[xm,ym,dyx]=linhist(gal.Z_SPEC,30);
line(xm,dyx/sum(ym),'color','g','linestyle','--');
xlim([0,0.6]);
ylabel('$dP/dz_{GAMA}$');
l=legend([h1,h2],'matched','all');set(l,'box','off','color','none');
% print('-depsc','GAMA_SRC_matchZ.eps');

x=0:0.02:1;
myfigure;
linhist(zsource(:,1),x,'stairs','r-');
hold on;
linhist([e{1}(:,9);e{2}(:,9);e{3}(:,9)],x,'stairs','r--');
plot(x+0.01,source_zdistr(x)*(x(2)-x(1))*numel([e{1}(:,9);e{2}(:,9);e{3}(:,9)]),'r-.');
linhist(gal.Z_SPEC(f),x,'stairs','g-');
linhist(gal.Z_SPEC,x,'stairs','g--');
l=legend('Matched','All','M05','GAMA Mactched','GAMA All');
set(l,'box','off');
xlabel('z');
ylabel('dN');
xlim([0,1]);
% print('-depsc','GAMA_SRC_matchZ_distr.eps');
%%
load gama_redshifts.mat
x=0:0.02:1;
myfigure;
[xm,ym,dyx]=linhist(zsdss(:,1),x);%,'stairsnorm','r-');
plot(xm,dyx/sum(ym),'r-');
hold on;
[xm,ym,dyx]=linhist(zsdss(:,3),x);%,'stairsnorm','r-');
plot(xm,dyx/sum(ym),'g-');
[xm,ym,dyx]=linhist(zsdss(:,5),x);%,'stairsnorm','r-');
plot(xm,dyx/sum(ym),'b-');
[xm,ym,dyx]=linhist(zsource(:,1),x);%,'stairsnorm','r-');
plot(xm,dyx/sum(ym),'y-');
[xm,ym,dyx]=linhist(zph4,x);%,'stairsnorm','r-');
plot(xm,dyx/sum(ym),'c-');
plot(x+0.01,source_zdistr(x),'m:');
[xm,ym,dyx]=linhist(zsp4,x);%,'stairsnorm','g-');
plot(xm,dyx/sum(ym),'k-','linewidth',2);
[xm,ym,dyx]=linhist(ztmpl,x);%,'stairsnorm','r--');
plot(xm,dyx/sum(ym),'r--');
[xm,ym,dyx]=linhist(zcc2,x);%,'stairsnorm','r--');
plot(xm,dyx/sum(ym),'g--');
[xm,ym,dyx]=linhist(zd1,x);%,'stairsnorm','r--');
plot(xm,dyx/sum(ym),'b--');
[xm,ym,dyx]=linhist([e{1}(:,9);e{2}(:,9);e{3}(:,9)],x);%,'stairs','r--');
plot(xm,dyx/sum(ym),'y--');
[xm,ym,dyx]=linhist(zsp4,x);%,'stairsnorm','g--');
plot(xm,dyx/sum(ym),'k--');
l=legend('sdss tmpl','sdss cc2','sdss d1','ZEBRA','gama photoz','M05','GAMA spec','All''s');
set(l,'box','off');
xlabel('z');
ylabel('dP/dz');
xlim([0,1]);
print('-depsc',['/home/kam/Projects/Lensing/data/srccat-gama/zgama_zsdss_distr_all.eps']);
%% magnitudes dependence of photoz
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
%% uppper and lower z-photo CI
myfigure;
plot(gal.Z_SPEC(f),zsource(:,4),'b.','markersize',3);hold on;
plot([0,0.6],[0,0.6],'r-');
axis([0,0.6,0,1]);
myfigure;
plot(gal.Z_SPEC(f),zsource(:,5),'g.','markersize',3);hold on;
plot([0,0.6],[0,0.6],'r-');
axis([0,0.6,0,1]);
%% error distribution
x=-2:0.02:2;
zerr=-[zsource(:,1)-zsource(:,2);zsource(:,1)-zsource(:,3)];
zerr2=-[zsource(:,1)-zsource(:,4);zsource(:,1)-zsource(:,5)];
myfigure;
[xm,ym,dyx]=linhist(zsource(:,1)-gal.Z_SPEC(f),x,'stairsnorm','r');
hold on;
[xmz,ymz,dyxz]=linhist(zerr,x,'stairsnorm','b');
[xmz2,ymz2,dyxz2]=linhist(zerr2,x,'stairsnorm','k');
% axis([0,0.6,0,1]);
% myfigure;plot(xm,dyx/sum(ym),'rx-',xmz,dyxz/sum(ymz),'ko-');
myfigure;stairs(x(1:end-1),dyx/sum(ym),'r');hold on;stairs(x(1:end-1),dyxz/sum(ymz),'b');stairs(x(1:end-1),dyxz2/sum(ymz2),'k');
[mu,sigma]=normfit(zsource(:,1)-gal.Z_SPEC(f));
[muz,sigmaz]=normfit(zerr);
[muz2,sigmaz2]=normfit(zerr2);
hold on;
plot(xm,normpdf(xm,mu,sigma),'r--');
plot(xmz,normpdf(xmz,muz,sigmaz),'b--');
plot(xmz2,normpdf(xmz2,muz2,sigmaz2),'k--');
xlim([-0.5,0.5]);
xlabel('$\Delta z$');ylabel('$dP/d\Delta z$');
l=legend('$z_{ph}-z_{spec}$','$z_{\sigma}-z_{ph}$','$z_{2\sigma}-z_{ph}$','Gaussian fits');
set(l,'box','off','interpreter','latex');
print('-depsc','srccat-gama/zerr_dist.eps');
%% template decomposition and magnitude comparison
myfigure;
linhist(ztemp(zsource(:,1)>0.36),-0.5:1:20.5,'stairsnorm','r');
hold on;
linhist(ztemp(zsource(:,1)<=0.36),-0.5:1:20.5,'stairsnorm','g');
l=legend('$z_{ph}>0.36$','$z_{ph}<0.36$');
set(l,'interpreter','latex');
xlabel('template');ylabel('Fraction');
% print('-depsc','photoz_templates.eps');

myfigure;
plot(gal.PETROMAG_R_SDSS(f),rsource,'.');
% axis([10,25,10,25]);
hold on;
plot([12,20],[12,20],'r-');
xlabel('$r_{pet}$');ylabel('$r_{mod}$');
% print('-depsc','GAMA_SRC_matchR.eps');

ematch=[e{1}(idmatch(idmatch(:,2)==1,1),:);e{2}(idmatch(idmatch(:,2)==2,1),:);e{3}(idmatch(idmatch(:,2)==3,1),:)];
%% template decomposition
xpos=0.07:0.3:1;
ypos=0.07:0.13:1;
ypos=ypos(end-1:-1:1);

zspec=gal.Z_SPEC(f);
h=zeros(7,3);
figure;
for i=0:20
    [x,y]=ind2sub([7,3],i+1);
    h(x,y)=axes('position',[xpos(y),ypos(x),0.3,0.13]);
ff=(ztemp==i);
plot(zspec(ff),zsource(ff,1),'.','markersize',3);
hold on;
axis([0,0.6,0,0.8]);
plot([0,0.6],[0,0.6],'r-');
l=legend(num2str(i));
set(l,'location','northwest','box','off');
end
set(h(:),'xtick',0:0.1:0.6,'xticklabel',0:0.1:0.5,'ytick',0:0.2:0.8,'yticklabel',0:0.2:0.6);
set(h(:,2),'yticklabel','');
set(h(:,3),'yaxislocation','right');
set(h(1,:),'xaxislocation','top');
set(h(2:6,:),'xticklabel','');
xlabel(h(7,2),'z_{GAMA}');
ylabel(h(4,1),'z_{ZEBRA}');
% title(num2str(i));
print('-depsc','GAMA_SRC_matchZ_temps.eps');