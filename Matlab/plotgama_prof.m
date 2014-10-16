function [h,r,s,es,coveff]=plotgama_prof(name,linespec,photosig)
% plot gama lensing profile with background subtraction and boost
% with subtraction and boost and photo-z error included
%   name: the mass bin to be ploted
%   linespec: line type
%   photosig: a flag controlling whether to use average photo-z signal or the original signal
%                    0: original
%                    1: photo-z averaged
%                    2: plot both

file=['/mnt/charon/Lensing/output_old/shearprof_',name,'.txt'];
data=importdata(file,'\t',1);
file=['/mnt/charon/Lensing/output_old/shearcov_',name,'.txt'];
cov=load(file);
coveff=cov;
for i=1:size(cov,1)
    for j=1:size(cov,2)
        coveff(i,j)=cov(i,j)/sqrt(cov(i,i))/sqrt(cov(j,j));
    end
end
for i=1:numel(data.textdata)
    if strcmp('r',data.textdata{i})
        rid=i;
    end
    if strcmp('p',data.textdata{i})
        pid=i;
    end
%     if strcmp('es',data.textdata{i})
%         esid=i;
%     end
    if strcmp('s',data.textdata{i})
        sid=i;
    end
    if strcmp('w',data.textdata{i})
        wid=i;
    end
end
file=['/mnt/charon/Lensing/output_old/randprof_',name,'.txt'];
rand=importdata(file,'\t',1);
for i=1:numel(rand.textdata)
     if strcmp('s',rand.textdata{i})
        randsid=i;
     end
     if strcmp('es',rand.textdata{i})
        randesid=i;
    end
    if strcmp('w',rand.textdata{i})
        randwid=i;
    end
    if strcmp('ew',rand.textdata{i})
        randewid=i;
    end
end
file=['/mnt/charon/Lensing/output_old/sphotoprof_',name,'.txt'];
sphoto=importdata(file,'\t',1);
for i=1:numel(sphoto.textdata)
    if strcmp('s',sphoto.textdata{i})
        sphotosid=i;
    end
    if strcmp('es',sphoto.textdata{i})
        sphotoesid=i;
    end
    if strcmp('w',sphoto.textdata{i})
        sphotowid=i;
    end
    if strcmp('ew',sphoto.textdata{i})
        sphotoewid=i;
    end
end
if nargin<2
    linespec='.';
end
if nargin<3
    photosig=0;
end

r=data.data(:,rid);

if photosig==0||photosig==2
p=data.data(:,pid);
s=data.data(:,sid)./p;     % QUESTION: which signal is better, real s or average s from random photo-z?
es=sqrt(diag(cov))./p;
%% photo-z error
es=es.^2+(sphoto.data(:,sphotoesid)./sphoto.data(:,1)).^2; % + photoz error
es=sqrt(es);
%% background subtraction
es=sqrt(es.^2+(rand.data(:,randesid)./rand.data(:,1)/sqrt(max(rand.data(:,end)))).^2);  % shear error + background error
s=s-rand.data(:,randsid)./rand.data(:,1);  %background subtraction
%% boost
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

es=sqrt(((es./s).^2+(ew./w).^2+(ewrand./wrand).^2).*(s.*w./wrand).^2); % + boost error
s=s.*w./wrand; %boost
%%
h=errorbar(r,signlog(s),signlog(s)-signlog(s-es),signlog(s+es)-signlog(s),linespec);
end

if photosig==1||photosig==2
%% ================using sphoto's s ===================
p=sphoto.data(:,1);
s=sphoto.data(:,sphotosid)./p;     % QUESTION: which signal is better, real s or average s from random photo-z?
es=sqrt(diag(cov))./data.data(:,pid);
%% photo-z error
es=es.^2+(sphoto.data(:,sphotoesid)./sphoto.data(:,1)).^2; % + photoz error
es=sqrt(es);
%% background subtraction
es=sqrt(es.^2+(rand.data(:,randesid)./rand.data(:,1)/sqrt(max(rand.data(:,end)))).^2);  % shear error + background error
s=s-rand.data(:,randsid)./rand.data(:,1);  %background subtraction
%% boost
% correct for averaging-bias
[w,ew]=weight_average_correction(sphoto.data(:,sphotowid),sphoto.data(:,sphotoewid),sphoto.data(:,end));
ew=ew/sqrt(max(sphoto.data(:,end)));  %this is not zero because it's not the original w used
% correct for averaging-bias
[wrand,ewrand]=weight_average_correction(rand.data(:,randwid),rand.data(:,randewid),rand.data(:,end));
ewrand=ewrand/sqrt(max(sphoto.data(:,end))); % this is also error on the average

es=sqrt(((es./s).^2+(ew./w).^2+(ewrand./wrand).^2).*(s.*w./wrand).^2); % + boost error
s=s.*w./wrand; %boost
%% 
hold on;
h=errorbar(r,signlog(s),signlog(s)-signlog(s-es),signlog(s+es)-signlog(s),'o');
end

set(gca,'xscale','log');
xlabel('r/(Mpc/h)');ylabel('$slog(\Delta\Sigma)$','interpreter','latex');