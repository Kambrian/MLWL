function load_GAMACore()
% to load GAMACoreDR1 data; if cols is specified, only load columns given
% by cols; otherwise load default IDs and Z

datadir='/home/kam/Projects/Lensing/data';
file=fullfile(datadir,'GamaCoreDR1_v1.csv');
fid=fopen(file);

%read header
txt=fgetl(fid);
head=textscan(txt,'%q','delimiter',',');
head=head{1};
ncol=numel(head);

%       Iid Gid Sid  ra dec r z  zQ  zS  zD  zsn zID p 
format='%s %d32 %u64 %f %f %f %f %d8 %d8 %d32 %f %s %s';
for i=14:33  %various magnitudes
    format=[format,' %f']; %#ok<AGROW>
end

%load data
data=textscan(fid,format,'delimiter',','); %#ok<NASGU>
fclose(fid);

%create variables
for i=1:ncol
    if ~isvarname(head{i})
        head{i}=genvarname(head{i});
    end
%     eval(['gama.',head{i},'=data{i};']);
    assignin('caller',head{i},data{i});
end

