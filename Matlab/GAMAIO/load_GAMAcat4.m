function gama=load_GAMAcat4(cols,formats)
%gama galaxy catalogue v4
datadir='/home/kam/Projects/Lensing/data';
file=fullfile(datadir,'catgama_v4.csv');
fid=fopen(file);

%read header
txt=fgetl(fid);
head=textscan(txt,'%q','delimiter',',');
head=head{1};
ncol=numel(head);
if ~issorted(cols)||cols(1)>cols(end)
    error('fields must be sorted ascending');
end

format='';
j=1;
for i=1:ncol  %various magnitudes
    if i==cols(j)
        format=[format,formats{j}];
        if j<numel(cols)
        j=j+1;
        end
    else
    format=[format,' %*s']; %#ok<AGROW>
    end
end

%load data
data=textscan(fid,format,1,'delimiter',','); %#ok<NASGU>
fclose(fid);

%create variables
for i=1:numel(cols)
    j=cols(i);
    if ~isvarname(head{j})
        head{j}=genvarname(head{j});
    end
    if nargout
    eval(['gama.',head{j},'=data{i};']);
    else
    assignin('caller',head{j},data{i});
    end
end

