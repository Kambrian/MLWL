% geometrical plots
%% group pie
figure;polar(cat09.CenRA./180*pi,cat09.MedianZ,'.r');
hold on;
polar(cat12.CenRA./180*pi,cat12.MedianZ,'.g');
polar(cat15.CenRA./180*pi,cat15.MedianZ,'.k');

polar(catmock.CenRA./180*pi+pi,catmock.MedianZ,'.b');
set(findobj('type','line'),'markersize',1);

text(0.7*cos(145/180*pi),0.7*sin(145/180*pi),'G09','color','r');
text(0.7*cos(180/180*pi),0.7*sin(180/180*pi),'G12','color','g');
text(0.7*cos(210/180*pi),0.7*sin(210/180*pi),'G015','color','k');
text(0.6*cos(0/180*pi),0.6*sin(0/180*pi),'MOCK','color','b');
title('GAMA Groups');
% 
% f=catmock.Mult>5;
% polar(catmock.CenRA(f),catmock.MedianZ(f),'.');
print('-depsc','GAMAgroup.eps');
%% galaxy pie
figure;polar(g09.RA./180*pi,g09.Z_SPEC,'.r');
hold on;
polar(g12.RA./180*pi,g12.Z_SPEC,'.g');
polar(g15.RA./180*pi,g15.Z_SPEC,'.k');

polar(mock.ra./180*pi+pi,mock.zz,'.b');
set(findobj('type','line'),'markersize',1);

text(0.7*cos(145/180*pi),0.7*sin(145/180*pi),'G09','color','r');
text(0.7*cos(180/180*pi),0.7*sin(180/180*pi),'G12','color','g');
text(0.7*cos(210/180*pi),0.7*sin(210/180*pi),'G015','color','k');
text(0.6*cos(0/180*pi),0.6*sin(0/180*pi),'MOCK','color','b');
title('GAMA Galaxies');

print('-depsc','GAMAgal.eps');
%% mock pie (redshift distortion)
figure;
polar(mock.ra./180*pi,mock.zz,'.r');
hold on;

polar(mock.ra./180*pi+pi,mock.zzhub,'.b');
set(findobj('type','line'),'markersize',1);

text(0.8*cos(180/180*pi),0.8*sin(180/180*pi),'Observed','color','r');
text(0.7*cos(0/180*pi),0.7*sin(0/180*pi),'Real','color','b');
title('Mock Galaxies');

print('-depsc','GAMAmock.eps');
%% colorful pie to show size? junk
colors=colormap(hot(max(catmock.Mult)));
x=colormap;
y=x(catmock.Mult/max(catmock.Mult)*64,:);
figure;colorpolar(cat09.CenRA,cat09.MedianZ,colors(cat09.Mult,:));
hold on;
polar(cat12.CenRA./180*pi,cat12.MedianZ,'.','markerfacecolor',colors(cat12.Mult))
polar(cat15.CenRA./180*pi,cat15.MedianZ,'.','markerfacecolor',colors(cat15.Mult))

f=catmock.Mult>5;
colorpolar(catmock.CenRA(f),catmock.MedianZ(f),y(f,:));
figure;polar(catmock.CenRA./180*pi,catmock.MedianZ,'.','markerfacecolor',colors(catmock.Mult))
%% test galID: the galID for BCGRef is the row number in galaxy file
i=find(mock.ra==catmock.BCGRA(1)&mock.dec==catmock.BCGDEC(1)&mock.amag_g==catmock.BCGmag(1));
disp([i,mock.idsel(i),mock.idgroup(i),mock.shortID(i)])
disp('galID='),disp(catmock.BCGRef(1))
j=find(mock.ra==catmock.BCGRA(2)&mock.dec==catmock.BCGDEC(2)&mock.amag_g==catmock.BCGmag(2));
disp([j,mock.idsel(j),mock.idgroup(j),mock.shortID(j)])
disp('galID='),disp(catmock.BCGRef(2))
%% mock projection
myfigure;
linktmp=linkmock;
galtmp=mock1;
cattmp=catmock;

grpid=1;
gids=linktmp(1:61,1);
plot(galtmp.ra(gids),galtmp.dec(gids),'.')
hold on;
plot(cattmp.BCGRA(1),cattmp.BCGDEC(1),'p','markersize',18);
% plot(cattmp.CenRA(1),cattmp.CenDEC(1),'o','markersize',18);
l2=plot_circle([cattmp.CenRA(grpid),cattmp.CenDEC(grpid)],cattmp.Rad100(grpid)/(cattmp.MedianDist(grpid)/(1+cattmp.MedianZ(grpid)))/pi*180);
% 
% grpid=2;
% gids=linktmp(61+(1:74),1);
% plot(galtmp.ra(gids),galtmp.dec(gids),'r.')
% hold on;
% plot(cattmp.BCGRA(2),cattmp.BCGDEC(2),'rp','markersize',5);
% % plot(cattmp.CenRA(2),cattmp.CenDEC(2),'ro','markersize',18);
% l2=plot_circle([cattmp.CenRA(grpid),cattmp.CenDEC(grpid)],cattmp.Rad100(grpid)/(cattmp.MedianDist(grpid)/(1+cattmp.MedianZ(grpid)))/pi*180);
% 
% grpid=3;
% gids=linktmp(sum(cattmp.Mult(1:2))+(1:cattmp.Mult(3)),1);
% plot(galtmp.ra(gids),galtmp.dec(gids),'g.')
% hold on;
% plot(cattmp.BCGRA(3),cattmp.BCGDEC(3),'gp','markersize',18);
% % plot(cattmp.CenRA(3),cattmp.CenDEC(3),'go','markersize',18);
% l2=plot_circle([cattmp.CenRA(grpid),cattmp.CenDEC(grpid)],cattmp.Rad100(grpid)/(cattmp.MedianDist(grpid)/(1+cattmp.MedianZ(grpid)))/pi*180);
% xlabel('RA');ylabel('DEC');
% print('-depsc','mockgroup3.eps');
axis equal
%% group projection
linktmp=link15;
galtmp=g15;
cattmp=cat15;
linkgrpid=[cattmp.Mult,(1:numel(cattmp.Mult))'];
linkgrpid=sortrows(linkgrpid,-1);%most abundant groups first

grpmax=500;
%map redshift to colors, zmin->zmax: blue->red;
colors=colormap(jet);
cstep=(max(cattmp.MedianZ(linkgrpid(1:grpmax,2)))-min(cattmp.MedianZ(linkgrpid(1:grpmax,2))))/64;
cmap=floor((cattmp.MedianZ-min(cattmp.MedianZ(linkgrpid(1:grpmax,2))))/cstep);
cmap(cmap==0)=1;cmap(cmap>64)=64;
% colors=colormap(jet(max(cattmp.Mult)));
% colors=colors(end:-1:1,:);

myfigure;
for i=1:grpmax
grpid=linkgrpid(i,2);
gids=linktmp(sum(cattmp.Mult(1:grpid-1))+(1:cattmp.Mult(grpid)),1);
for j=1:numel(gids), gids(j)=find(galtmp.CATA_INDEX==gids(j));end
l1=plot(galtmp.RA(gids),galtmp.DEC(gids),'.','markersize',3);
hold on;
l2=plot_circle([cattmp.CenRA(grpid),cattmp.CenDEC(grpid)],cattmp.Rad100(grpid)/(cattmp.MedianDist(grpid)/(1+cattmp.MedianZ(grpid)))/pi*180);
% l2=plot(cattmp.BCGRA(grpid),cattmp.BCGDEC(grpid),'o','markersize',cattmp.Mult(grpid)^(1/3)*5);
% set([l1,l2],'color',colors(linkgrpid(i,1),:));
set([l1,l2],'color',colors(cmap(linkgrpid(i,2)),:));
end
axis equal
% colorbar;

print('-depsc','G15groups.eps');