function out=loadGAMAcsv(file,version,nskip)
% myfile to load a csv file for GAMA data
% inputs: 
%           file: filename
%           version: GAMA file version; default to 2 if not specified
%outputs:
%           if not specified, output all columns into calling function
%           otherwise return a struct to out

% file='mockcutcatvol1.csv';
fid=fopen(file);
if nargin<2
    version=2;
end
if nargin<3
    nskip=0;
end

if 2==version
delim=',';
headformat='%q';
LongIDs='idgroup';
end

if 4==version
    delim='';
    headformat='%s';
    LongIDs='MillSimID';
end

if 6==version
    delim=',';
    headformat='%s';
    LongIDs={'UniqueID','galID','idhalo','IterRef','BCGRef'};
end

%read header
for i=1:nskip
    fgetl(fid);
end

txt=fgetl(fid);
if isempty(delim)
    head=textscan(txt,headformat);
else
    head=textscan(txt,headformat,'delimiter',delim);
end
head=head{1}; %open the cell; all the members are still cells

%construct file format specifier
ncol=numel(head);
format=[];
for i=1:ncol
    if is_in(head{i},LongIDs)
%         txtcol=i;
        if 2==version
            format=[format,sprintf('\"'),'%u64',sprintf('\"'),' '];
        elseif 6==version
            if strcmp(head{i},'UniqueID')
                format=[format,'%s',' '];
            else
                format=[format,'%u64', ' '];
            end
        else
             format=[format,'%u64',' '];
        end
    else
        format=[format,'%f ']; %#ok<*AGROW>
    end
end
format(end)=[];

%load data
if isempty(delim)
    data=textscan(fid,format,'TreatAsEmpty','NA'); %#ok<NASGU>
else
    data=textscan(fid,format,'delimiter',delim,'TreatAsEmpty','NA'); %#ok<NASGU>
end
fclose(fid);

%create variables
for i=1:ncol
    if ~isvarname(head{i})
        head{i}=genvarname(head{i});
    end
    eval(['out.',head{i},'=data{i};']);
% assignin('caller',head{i},data{i});
end

    function y=is_in(el,set)
        y=0;
        for x=set
            if strcmp(el,x)
                y=1;
                return;
            end
        end
    end
end
