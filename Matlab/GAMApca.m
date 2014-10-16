 cat=structcat([cat{1};cat{2};cat{3}]);
names=fieldnames(cat);
catmat=zeros(numel(cat.(names{1})),numel(names));
for i=1:numel(names)
    catmat(:,i)=cat.(names{i});
end

catmat=[cat.Mult,cat.MedianZ,cat.Rad1Sig,cat.VelDisp,cat.Mass,cat.LumMass,cat.BCGmag,cat.TotFluxInt,cat.TotSMInt];
mmat=[cat.Mass,cat.VelDisp,cat.Rad1Sig];

for i=1:3
mmat(mmat(:,i)<=0,:)=[];
end

lmat=log(mmat);
