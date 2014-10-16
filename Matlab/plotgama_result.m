function [h,r,y,ey]=plotgama_result(name,type,linespec)
%plot s or n profile
% this s is raw data, only corrected for responsivity

file=['/mnt/charon/Lensing/output/shearprof_',name,'.txt'];
data=importdata(file,'\t',1);
file=['/mnt/charon/Lensing/output/shearcov_',name,'.txt'];
cov=load(file);
for i=1:numel(data.textdata)
    if strcmp('r',data.textdata{i})
        rid=i;
    end
    if strcmp('p',data.textdata{i})
        pid=i;
    end
    if strcmp('es',data.textdata{i})
        esid=i;
    end
    if strcmp('s',data.textdata{i})
        sid=i;
    end
     if strcmp('n',data.textdata{i})
        nid=i;
    end
end

if nargin<3
    linespec='.';
end
switch type
    case 's'
p=data.data(:,pid);        
s=data.data(:,sid)./p;
% es=data.data(:,esid)./p;
es=sqrt(diag(cov))./p;
r=data.data(:,rid);y=s;ey=es;
h=errorbar(data.data(:,rid),signlog(s),signlog(s)-signlog(s-es),signlog(s+es)-signlog(s),linespec);
set(gca,'xscale','log');
% set(gca,'yscale','log');
    case 'n'
    file=['/mnt/charon/Lensing/output/randprof_',name,'.txt'];
    rand=importdata(file,'\t',1);    
    for i=1:numel(rand.textdata)
    if strcmp('n',rand.textdata{i})
        randnid=i;
    end
    if strcmp('en',rand.textdata{i})
        randenid=i;
    end
    end
    [randn,randen]=weight_average_correction(rand.data(:,randnid),rand.data(:,randenid),rand.data(:,end));%single measurement error
    n=data.data(:,nid);
    en=0.;
    en=sqrt((en./n).^2+(randen./randn).^2).*abs(n./randn);
    n=n./randn;
% en=sqrt((1./data.data(:,nid)+(rand.data(:,randenid)./sqrt(rand.data(:,end))./rand.data(:,randnid)).^2)).*abs(n);
%=======
%  rand.data(:,8)./sqrt(rand.data(:,4))
%
%  error in n can be approximated to be poisson on small scales, when the
%  src galaxy counts do not involve duplicate entries; on large scales,
%  poisson error far underestimates the real scatter
%=========
h=errorbar(data.data(:,rid),n,en,linespec);
r=data.data(:,rid);y=n;ey=en;
set(gca,'xscale','log');
end
    