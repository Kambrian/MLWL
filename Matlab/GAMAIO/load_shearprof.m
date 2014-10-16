function data=load_shearprof(file,mbin,zbin,rbin)

data=struct('s',{},'w',{},'es',{},'r',{},'r2',{},'n',{},'nstack',{});
fid=fopen(file);
tmp=fread(fid,1,'int32');
for i=1:mbin
    for j=1:zbin
        data(i,j).s=fread(fid,rbin,'float64');
        data(i,j).w=fread(fid,rbin,'float64');
        data(i,j).es=fread(fid,rbin,'float64');
        data(i,j).r=fread(fid,rbin,'float64');
        data(i,j).r2=fread(fid,rbin,'float64');
        data(i,j).n=fread(fid,rbin,'int32');
        data(i,j).nstack=fread(fid,1,'int32');
    end
end
tmp2=fread(fid,1,'int32');
if(tmp~=tmp2)
    error('file corruption or wrong');
end
fclose(fid);