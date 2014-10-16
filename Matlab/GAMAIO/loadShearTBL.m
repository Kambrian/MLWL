function data=loadShearTBL(n,ind)
% to load shear map for gama provided by Mandelbaum
% n: gama region id; 
% ind: columns to load
datadir='/home/kam/Projects/Lensing/data/srccat-gama';
file=fullfile(datadir,['srccat-010-gama-',num2str(n),'.tbl']);
fp=fopen(file);
      %[1   2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 12 22 23 24 25 26 27 28 29 30 ]
format=['%f %f %f %d %d %d %d %d %d %d %s %s %s %f %f %f %f %f %f %f %f %f %f %d %d %d %f %f %f %f ',...
        '%f %d %d %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f %f %f %f'];
tbl=textscan(fp,format);
fclose(fp);
tbl{11}=hex2dec(tbl{11});
tbl{12}=hex2dec(tbl{12});
tbl{13}=hex2dec(tbl{13});
if nargin<2, ind=1:60; end
data=zeros(numel(tbl{1}),numel(ind));
for i=1:numel(ind)
data(:,i)=tbl{ind(i)};
end
% data=cell2mat(tbl);