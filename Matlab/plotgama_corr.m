function [h,r,w,ew]=plotgama_corr(name,linespec)
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
    if strcmp('w',data.textdata{i})
        wid=i;
    end
end
file=['/mnt/charon/Lensing/output/randprof_',name,'.txt'];
rand=importdata(file,'\t',1);
for i=1:numel(rand.textdata)
     if strcmp('s',rand.textdata{i})
        randsid=i;
    end
    if strcmp('w',rand.textdata{i})
        randwid=i;
    end
    if strcmp('ew',rand.textdata{i})
        randewid=i;
    end
end


p=data.data(:,pid);
r=data.data(:,rid);

w=data.data(:,wid);
ew=0.;  % this should be set zero because w is cancled out from
                 % sum(w*e)/sum(w)*sum(w)/sum(w_rand), so each w is
                 % effectively a constant
% correct for averaging-bias
[wrand,ewrand]=weight_average_correction(rand.data(:,randwid),rand.data(:,randewid),rand.data(:,end));
% here ewrand should be error on single measurement of the background, not
% error on the mean; because the correction is expected to be done using
% w/wrand with w and wrand measured once at the same time. we do not have
% wrand measured together with w, so we use the average wrand, and quote
% the uncertaity for a single measurement 

ew=sqrt((ew./w).^2+(ewrand./wrand).^2).*abs(w./wrand);
w=w./wrand;
% es./s
% ew./w
% s=s-rand.data(:,randsid);
% es=sqrt((es./s).^2+(ew./w).^2).*abs(s.*w);
% s=s.*w;

if nargin<2
    linespec='.';
end

h=errorbar(r,w,ew,linespec);

set(gca,'xscale','log');
