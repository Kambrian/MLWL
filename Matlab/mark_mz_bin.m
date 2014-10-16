function labels=mark_mz_bin(bintype)
% mark bin ids onto groups

global grp

labels=cell(4,1);

for sky=1:3
grp{sky}.zid=zeros(numel(grp{sky}.GroupID),1);
grp{sky}.zid(grp{sky}.Zfof<=0.2)=1;
grp{sky}.zid(grp{sky}.Zfof>0.2)=2;

grp{sky}.mid=zeros(numel(grp{sky}.GroupID),1);
switch bintype
    case 'halomass'
mcut=[1e12,1e13,1e14];
grp{sky}.mid(grp{sky}.HybMass<=mcut(1))=1;
grp{sky}.mid(grp{sky}.HybMass>mcut(1)&grp{sky}.HybMass<=mcut(2))=2;
grp{sky}.mid(grp{sky}.HybMass>mcut(2)&grp{sky}.HybMass<=mcut(3))=3;
grp{sky}.mid(grp{sky}.HybMass>mcut(3))=4;
labels{1}=['$M_h\leq',printexp10(mcut(1)),'M_\odot/h$'];
labels{2}=['$',printexp10(mcut(1)),'M_\odot/h<M_h\leq',printexp10(mcut(2)),'M_\odot/h$'];
labels{3}=['$',printexp10(mcut(2)),'M_\odot/h<M_h\leq',printexp10(mcut(3)),'M_\odot/h$'];
labels{4}=['$M_h>',printexp10(mcut(3)),'M_\odot/h$'];
    case 'lummass'
mcut=[5e11,1e12,2e12];
grp{sky}.mid(grp{sky}.LumMass<=mcut(1))=1;
grp{sky}.mid(grp{sky}.LumMass>mcut(1)&grp{sky}.LumMass<=mcut(2))=2;
grp{sky}.mid(grp{sky}.LumMass>mcut(2)&grp{sky}.LumMass<=mcut(3))=3;
grp{sky}.mid(grp{sky}.LumMass>mcut(3))=4;
labels{1}=['$M_L\leq',printexp10(mcut(1)),'M_\odot/h$'];
labels{2}=['$',printexp10(mcut(1)),'M_\odot/h<M_L\leq',printexp10(mcut(2)),'M_\odot/h$'];
labels{3}=['$',printexp10(mcut(2)),'M_\odot/h<M_L\leq',printexp10(mcut(3)),'M_\odot/h$'];
labels{4}=['$M_L>',printexp10(mcut(3)),'M_\odot/h$'];
    case 'dynmass'
mcut=[2e12,2e13,1e14];
grp{sky}.mid(grp{sky}.MassProxy<=mcut(1))=1;
grp{sky}.mid(grp{sky}.MassProxy>mcut(1)&grp{sky}.MassProxy<=mcut(2))=2;
grp{sky}.mid(grp{sky}.MassProxy>mcut(2)&grp{sky}.MassProxy<=mcut(3))=3;
grp{sky}.mid(grp{sky}.MassProxy>mcut(3))=4;
labels{1}=['$M_d\leq',printexp10(mcut(1)),'M_\odot/h$'];
labels{2}=['$',printexp10(mcut(1)),'M_\odot/h<M_d\leq',printexp10(mcut(2)),'M_\odot/h$'];
labels{3}=['$',printexp10(mcut(2)),'M_\odot/h<M_d\leq',printexp10(mcut(3)),'M_\odot/h$'];
labels{4}=['$M_d>',printexp10(mcut(3)),'M_\odot/h$'];        
    case 'bcgmag'
mcut=[-20.5,-21,-21.5];
grp{sky}.mid(grp{sky}.BCGmag>=mcut(1))=1;
grp{sky}.mid(grp{sky}.BCGmag<mcut(1)&grp{sky}.BCGmag>=mcut(2))=2;
grp{sky}.mid(grp{sky}.BCGmag<mcut(2)&grp{sky}.BCGmag>=mcut(3))=3;
grp{sky}.mid(grp{sky}.BCGmag<mcut(3))=4;
labels{1}=['$M_{BCG}\geq',num2str(mcut(1)),'$'];
labels{2}=['$',num2str(mcut(1)),'>M_{BCG}\geq',num2str(mcut(2)),'$'];
labels{3}=['$',num2str(mcut(2)),'>M_{BCG}\geq',num2str(mcut(3)),'$'];
labels{4}=['$M_{BCG}<',num2str(mcut(3)),'$'];           
    case 'comvsize'
mcut=[0.05,0.1,0.2];
grp{sky}.mid(grp{sky}.Rad1Sig<=mcut(1))=1;
grp{sky}.mid(grp{sky}.Rad1Sig>mcut(1)&grp{sky}.Rad1Sig<=mcut(2))=2;
grp{sky}.mid(grp{sky}.Rad1Sig>mcut(2)&grp{sky}.Rad1Sig<=mcut(3))=3;
grp{sky}.mid(grp{sky}.Rad1Sig>mcut(3))=4;
labels{1}=['$R\leq',printexp10(mcut(1)),'Mpc/h$'];
labels{2}=['$',printexp10(mcut(1)),'Mpc/h<R\leq',printexp10(mcut(2)),'Mpc/h$'];
labels{3}=['$',printexp10(mcut(2)),'Mpc/h<R\leq',printexp10(mcut(3)),'Mpc/h$'];
labels{4}=['$R>',printexp10(mcut(3)),'Mpc/h$'];    
    otherwise
        error('unknown bin type');
end
end