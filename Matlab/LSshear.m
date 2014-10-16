rmin=100;
cd /home/kam/Projects/Lensing/data
cat{1}=loadGAMAcsv('groupsv2/groupcatv2G09opt194lim194.csv');
cat{2}=loadGAMAcsv('groupsv2/groupcatv2G12opt194lim194.csv');
cat{3}=loadGAMAcsv('groupsv2/groupcatv2G15opt194lim194.csv');
catall=structcat([cat{1},cat{2},cat{3}]);
cd /home/kam/Projects/Lensing/output
ngrp=numel(catall.Mult);
fp=fopen(['Signal_LS',num2str(rmin),'.dat']);
signal=fread(fp,inf,'float32');
fclose(fp);
fp=fopen(['Noise_LS',num2str(rmin),'.dat']);
noise=fread(fp,inf,'float32');
fclose(fp);
if ngrp~=numel(signal)
    error('ngrp mismatch');
end
n=zeros(3,1);
s=cell(3,1);
j=0;
for i=1:3
    n(i)=numel(cat{i}.Mult);
    s{i}=signal(j+(1:n(i)));
    j=j+n(i);
end

%%
myfigure;
for i=1:3
    ra=cat{i}.BCGRA;
    dec=cat{i}.BCGDEC;
    xi=linspace(min(ra),max(ra),100);
    yi=linspace(min(dec),max(dec),100);
    F=TriScatteredInterp(ra,dec,signlog(s{i}));
    [xx,yy]=meshgrid(xi,yi);
    zz=F(xx,yy);
    %  [xx,yy,zz]=griddata(ra,dec,signlog(s{i}),xi,yi');
    subplot(3,1,i);
    imagesc(minmax(xi),minmax(yi),zz);
    set(gca,'ydir','normal');
    % contour(xx,yy,zz);
end
print('-depsc',['signal_imageLS',num2str(rmin),'.eps']);
%%
ncolor=256;
myfigure;
c=colormap(jet(ncolor));
for i=1:3
    ra=cat{i}.BCGRA;
    dec=cat{i}.BCGDEC;
    data=signlog(s{i});
    ra(isnan(data))=[];
    dec(isnan(data))=[];
    data(isnan(data))=[];
%     [tmp,ic]=histc(data,linspace(min(data)-eps,max(data)+eps,ncolor));
     [tmp,ic]=histc(data,linspace(-5,5,ncolor));    ic(data>5)=ncolor;ic(data<-5)=1;
    subplot(3,1,i);
    for j=1:numel(ra)
    plot(ra(j),dec(j),'.','color',c(ic(j),:));
    hold on;
    end
    axis tight;
end
print('-depsc',['signal_mapLS',num2str(rmin),'.eps']);
%%
myfigure;
for i=1:3
    ra=cat{i}.BCGRA;
    dec=cat{i}.BCGDEC;
    f=s{i}>0;
    subplot(3,1,i);
    plot(ra(f),dec(f),'r.');
    hold on;
    plot(ra(~f),dec(~f),'g.');
    axis tight;
end
% print('-depsc',['signalsign_mapLS',num2str(rmin),'.eps']);
% axis equal;
%%
figure;scatter3(catall.BCGRA,catall.BCGDEC,signlog(signal));
figure;
hist(abs(signal)./noise)