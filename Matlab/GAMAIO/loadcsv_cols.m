function datacat=loadcsv_cols(file,cols,formats)
% to load the csv file with one header line
% input: file: filename to be read
%        cols: vector, specifying the columns to be read
%        formats:cell array of the same size as cols, specifying the format
%                of each column in cols
% output: datacat (optional): a struct containing the columns loaded;
%         if output is not assigned, the data loaded would be directly
%         exported to the memory space of the calling function
%
% J.X. Han, 10/02/2010
% 

% datadir='/home/kam/Projects/Lensing/data';
% file=fullfile(datadir,'catgama_v4.csv');
fid=fopen(file);

%read header
txt=fgetl(fid);
head=textscan(txt,'%q','delimiter',',');
head=head{1};
ncol=numel(head);
[cols,IDX]=sort(cols);
formats=formats(IDX);

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
data=textscan(fid,format,'delimiter',','); %#ok<NASGU>
fclose(fid);

%create variables
for i=1:numel(cols)
    j=cols(i);
    if ~isvarname(head{j})
        head{j}=genvarname(head{j});
    end
    if nargout
    eval(['datacat.',head{j},'=data{i};']);
    else
    assignin('caller',head{j},data{i});
    end
end

