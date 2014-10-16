function a=shearprof_rebin(file,nbinvar)
data=load_shearprof(file,4,2,60);
% data=load_shearprof('shearprof_rand1.dat');
data2=struct([]);
for i=1:4
    for j=1:2
        data(i,j).s(isnan(data(i,j).s))=0;
        data(i,j).r(isnan(data(i,j).r))=0;
        data(i,j).r2(isnan(data(i,j).r2))=0;
        data(i,j).es(isnan(data(i,j).es))=0;
    end
end
for i=1:4
    data2(i).nstack=data(i,1).nstack+data(i,2).nstack;
    data2(i).s=data(i,1).s.*data(i,1).w+data(i,2).s.*data(i,2).w;
    data2(i).es=sqrt((data(i,1).es.*data(i,1).w).^2+(data(i,2).es.*data(i,2).w).^2);
    data2(i).r=data(i,1).r.*data(i,1).w+data(i,2).r.*data(i,2).w;
    data2(i).r2=data(i,1).r2.*data(i,1).w+data(i,2).r2.*data(i,2).w;
    data2(i).w=data(i,1).w+data(i,2).w;
    data2(i).n=data(i,1).n+data(i,2).n;
    data2(i).s=data2(i).s./data2(i).w;
    data2(i).es=data2(i).es./data2(i).w;
    data2(i).r=data2(i).r./data2(i).w;
    data2(i).r2=data2(i).r2./data2(i).w;
    data2(i).s(isnan(data2(i).s))=0;
    data2(i).r(isnan(data2(i).r))=0;
    data2(i).r2(isnan(data2(i).r2))=0;
    data2(i).es(isnan(data2(i).es))=0;
end

% nbinvar=[3,6,12,15];
nbin=60;
a.rvar=cell(4,1);
a.ervar=cell(4,1);
a.svar=cell(4,1);
a.esvar=cell(4,1);
a.nvar=cell(4,1);
a.wvar=cell(4,1);
a.ngrps=zeros(4,1);
for i=1:4
    a.wvar{i}=sum(reshape(data2(i).w,nbin/nbinvar(i),nbinvar(i)),1);
    a.rvar{i}=sum(reshape(data2(i).r.*data2(i).w,nbin/nbinvar(i),nbinvar(i)),1)./sum(reshape(data2(i).w,nbin/nbinvar(i),nbinvar(i)),1);
    a.ervar{i}=sum(reshape(data2(i).r2.*data2(i).w,nbin/nbinvar(i),nbinvar(i)),1)./sum(reshape(data2(i).w,nbin/nbinvar(i),nbinvar(i)),1);
    a.nvar{i}=sum(reshape(data2(i).n,nbin/nbinvar(i),nbinvar(i)),1);
    a.svar{i}=sum(reshape(data2(i).s.*data2(i).w,nbin/nbinvar(i),nbinvar(i)),1)./sum(reshape(data2(i).w,nbin/nbinvar(i),nbinvar(i)),1);
    a.esvar{i}=sqrt(sum(reshape((data2(i).es.*data2(i).w).^2,nbin/nbinvar(i),nbinvar(i)),1))./sum(reshape(data2(i).w,nbin/nbinvar(i),nbinvar(i)),1);
    a.rvar{i}(isnan(a.rvar{i}))=0;
    a.ervar{i}(isnan(a.ervar{i}))=0;
    a.svar{i}(isnan(a.svar{i}))=0;
    a.esvar{i}(isnan(a.esvar{i}))=0;
    a.ngrps(i)=data2(i).nstack;
end

%%
% rmin=0.02;rmax=300;
% a.v=cell(4,1);
% for i=1:4
%      tmp=logspace(log10(rmin),log10(rmax),nbinvar(i)+1);
%     a.v{i}=pi*diff(tmp.^2);
% end
% 
% a.nmean=zeros(1,4);
% a.y=cell(4,1);a.ey=cell(4,1);
% for i=1:4
%     a.nmean(i)=sum(a.nvar{i})./rmax^2/pi/a.ngrps(i);
% %     data.nmean(i)=data.nmean(i)/nz_factor(data.zmean(i))*(pi/180/60)^2;
%     a.y{i}=(a.nvar{i}./a.v{i})/a.ngrps(i);
%     a.ey{i}=sqrt(a.nvar{i})./a.v{i}/a.ngrps(i);
% end
% 
% ltyp={'ro-','ro--';'g>-','g>--';'bs-','bs--';'kd-','kd--';};
% ldtyp={'ro:';'g>:';'bs:';'kd:';};
% lltyp={'ro-';'g>-';'bs-';'kd-';};
% mtyp={'ro';'g>';'bs';'kd';};
% color={'r';'g';'b';'k'};
% dir='/home/kam/Projects/Lensing/output/';
% 
% myfigure;
% for i=1:4
%     y=a.y{i};
%     ey=a.ey{i};
% %     y=y/nz_factor(datamock.zmean(i))*(pi/180/60)^2;
% %     ey=ey/nz_factor(datamock.zmean(i))*(pi/180/60)^2;
%     h1=errorbar(a.rvar{i}',y,ey,ldtyp{i});hold on;
% %     y=data.y{i};
% %     ey=data.ey{i};
% %     y=y/nz_factor(data.zmean(i))*(pi/180/60)^2;
% %     ey=ey/nz_factor(data.zmean(i))*(pi/180/60)^2;
% %     h2=errorbar(data.rvar{i}',y,ey,lltyp{i},'markerfacecolor',color{i});hold on;
% end
% set(gca,'xscale','log','yscale','log');