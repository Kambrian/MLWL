function data=load_randprof(file,mbin,zbin,rbin)

data=struct('s',{},'w',{},'n',{},'es',{},'ew',{},'en',{},'nrand',{});
fid=fopen(file);
tmp=fread(fid,1,'int32');
for i=1:mbin
    for j=1:zbin
        data(i,j).s=fread(fid,rbin,'float64');
        data(i,j).w=fread(fid,rbin,'float64');
        data(i,j).n=fread(fid,rbin,'float64');
        data(i,j).es=fread(fid,rbin,'float64');
        data(i,j).ew=fread(fid,rbin,'float64');
        data(i,j).en=fread(fid,rbin,'float64');
        data(i,j).nrand=fread(fid,1,'int32');
    end
end
tmp2=fread(fid,1,'int32');
if(isempty(tmp2)||tmp~=tmp2)
    error('file corruption or wrong');
end
fclose(fid);