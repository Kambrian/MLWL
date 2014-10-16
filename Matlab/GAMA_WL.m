function data=GAMA_WL(macro)
%% switches:
% log/lin scale
% BCGRA\DEC, CENRA\DEC
% Multiplicity filter: N==2,3 include or not? N==2 further split into two populations
% Mass Proxy, bin division
% rmin, rmax
% physical/comoving binning
% error estimation method (error propagation; photo-z error inclusion;
%                          bootstrap;)
% photo-z accuracy filter
% estimators: M_zeta, M_zetac, DSig, Dsig_Inv, Mandelbaum

% to investigate:****
% filtering trade-off: boost or decrease S/N, how to maximize
%%
GAMA_WL_init
%%
% tdata=cell(nbinm,1);
flag_comoving=macros.STACK_COMOV;
for sky=1:3
    ngrps=numel(grp{sky}.LumMass);
    prog=0;
    fprintf(1,'sky %02d:    ',sky);
    for gid=1:ngrps
        if gid==109
            109;
        end
        if gid>ngrps/10*prog, fprintf(1,['\b\b\b',num2str(prog*10,'%02d'),'%%']);prog=prog+1; end
        if grp{sky}.Mult(gid)<macros.MULT_MIN, continue; end   %additional filter
        mid=grp{sky}.mid(gid);zid=grp{sky}.zid(gid);
        ngrps_used(mid,zid)=ngrps_used(mid,zid)+1;
        [rsum,ersum,m,n,s,w,em,es,s_pred]=group_WL_signal(grpcenters{sky}(gid,:),[rmin,rmax(mid,zid)],grp{sky}.Zfof(gid),sky,nbin,flag_comoving,grp{sky}.LumMass(gid));
        %[rsum,ersum,m,n,s,w,em,es]=group_WL_signal_angular(grpcenters{sky}(gid,:),[rmin,rmax(mid,zid)],grp{sky}.Zfof(gid),sky,nbin,flag_comoving);
%         [et,ex,es,sigmc]=group_WL_test([grp{sky}.CenRA(gid),grp{sky}.CenDEC(gid)],[rmin,rmax(mid,zid)],grp{sky}.MedianZ(gid),sky,flag_comoving);
%         tdata{mid}=[tdata{mid};et,ex,es,sigmc];
        rr(:,mid,zid)=rr(:,mid,zid)+rsum;
        rr2(:,mid,zid)=rr2(:,mid,zid)+ersum;
        mr(:,mid,zid)=mr(:,mid,zid)+m;
        nr(:,mid,zid)=nr(:,mid,zid)+n;
        sr(:,mid,zid)=sr(:,mid,zid)+s;
        wr(:,mid,zid)=wr(:,mid,zid)+w;
        emr(:,mid,zid)=emr(:,mid,zid)+em;
        esr(:,mid,zid)=esr(:,mid,zid)+es;
        sr_pred(:,mid,zid)=sr_pred(:,mid,zid)+s_pred;
    end
    fprintf(1,'\n');
end
  rr=rr./wr;
  rr2=rr2./wr;
  mr=mr./nr/2;  %factor 2 comes form converting from chi to epsilon
  sr=sr./wr/2;
  sr_pred=sr_pred./wr/2;
  emr=sqrt(emr)./nr/2;  % conversion as well
  esr=sqrt(esr)./wr/2;
  
  rl=mr;
  for i=1:nbinm
      for j=1:nbinz
%           tmp=linspace(rmin,rmax(i,j),nbin+1);
          tmp=logspace(log10(rmin),log10(rmax(i,j)),nbin+1);
          rl(:,i,j)=tmp(1:end-1);
      end
  end
  mr(isnan(mr))=0;sr(isnan(sr))=0;emr(isnan(emr))=0;esr(isnan(esr))=0;sr_pred(isnan(sr_pred))=0;
  rr(isnan(rr))=0;rr2(isnan(rr2))=0;
  % rebin the redshifts into one bin
  m=sum(mr.*nr,3)./sum(nr,3);
  s=sum(sr.*wr,3)./sum(wr,3);
  s_pred=sum(sr_pred.*wr,3)./sum(wr,3);
  em=sqrt(sum((emr.*nr).^2,3))./sum(nr,3);
  es=sqrt(sum((esr.*wr).^2,3))./sum(wr,3);
  r=sum(rr.*wr,3)./sum(wr,3);
  r2=sum(rr2.*wr,3)./sum(wr,3);
  m(isnan(m))=0;s(isnan(s))=0;em(isnan(em))=0;es(isnan(es))=0;
  r(isnan(r))=0;r2(isnan(r2))=0;
  disp('all done')
%%
disp('saving...')
eval(['save ',gamaWL_fname(macros),'.mat mr nr sr wr rl rr rr2 rmax r r2 m s emr esr em es ngrps_used sr_pred s_pred labels macros'])
% save angular_shear.mat  mr nr sr wr rl rr rr2 rmax r r2 m s emr esr em es ngrps_used labels macros
% msgbox('GAMA job done');
% data=load([gamaWL_fname(macros),'.mat']);
data=load(['angular_shear.mat']);