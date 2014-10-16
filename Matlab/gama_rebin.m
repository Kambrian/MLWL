function a=gama_rebin(macro,nbinvar)

load(['/home/kam/Projects/Lensing/output/',gamaWL_fname(macro),'.mat']);
nbin=size(nr,1);

a.rvar=cell(4,1);
a.ervar=cell(4,1);
a.mvar=cell(4,1);
a.svar=cell(4,1);
a.spredvar=cell(4,1);
a.emvar=cell(4,1);
a.esvar=cell(4,1);
a.nvar=cell(4,1);
a.wvar=cell(4,1);
a.ngrps=sum(ngrps_used,2);
n=sum(nr,3);w=sum(wr,3);
n=[n;0 0 0 0];
for i=1:4
    a.wvar{i}=sum(reshape(w(:,i),nbin/nbinvar(i),nbinvar(i)));
    a.rvar{i}=sum(reshape(r(:,i).*w(:,i),nbin/nbinvar(i),nbinvar(i)))./sum(reshape(w(:,i),nbin/nbinvar(i),nbinvar(i)));
    a.ervar{i}=sum(reshape(r2(:,i).*w(:,i),nbin/nbinvar(i),nbinvar(i)))./sum(reshape(w(:,i),nbin/nbinvar(i),nbinvar(i)));
    a.rlvar{i}=rl(1:nbin/nbinvar(i):nbin,i,1);
    a.mvar{i}=m(1:nbin/nbinvar(i):nbin,i);  % resample rather than smooth in radius
    a.nvar{i}=n(1:nbin/nbinvar(i):nbin,i)-n(1+nbin/nbinvar(i):nbin/nbinvar(i):nbin+1,i);
    a.svar{i}=sum(reshape(s(:,i).*w(:,i),nbin/nbinvar(i),nbinvar(i)))./sum(reshape(w(:,i),nbin/nbinvar(i),nbinvar(i)));
    a.spredvar{i}=sum(reshape(s_pred(:,i).*w(:,i),nbin/nbinvar(i),nbinvar(i)))./sum(reshape(w(:,i),nbin/nbinvar(i),nbinvar(i)));
    a.emvar{i}=em(1:nbin/nbinvar(i):nbin,i);
    a.esvar{i}=sqrt(sum(reshape((es(:,i).*w(:,i)).^2,nbin/nbinvar(i),nbinvar(i))))./sum(reshape(w(:,i),nbin/nbinvar(i),nbinvar(i)));
    a.rvar{i}(isnan(a.rvar{i}))=0;
    a.ervar{i}(isnan(a.ervar{i}))=0;
    a.svar{i}(isnan(a.svar{i}))=0;
    a.spredvar{i}(isnan(a.spredvar{i}))=0;
    a.esvar{i}(isnan(a.esvar{i}))=0;
    a.mvar{i}(isnan(a.mvar{i}))=0;
    a.emvar{i}(isnan(a.emvar{i}))=0;
end

