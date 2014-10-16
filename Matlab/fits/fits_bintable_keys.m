function [header,fileinfo]=fits_bintable_keys(file,tableid)
% function [header,fileinfo]=fits_bintable_keys(file,tableid)
%
% collect column names from binary table of a fits files
% can return fitsinfo as well
% 
% file: fits filename
% tableid: bintable tableid, in case there are more than one bintables
%            available (starting from 1)
%
% ==Jiaxin Han, 11/02/2010, Durham ==
% Revision: 06/16/2011, Durham
%

if nargin<2, tableid=1; end

a=fitsinfo(file);
Keys=a.BinaryTable(tableid).Keywords;
header=cell(1,a.BinaryTable(tableid).NFields);
j=0;
for i=1:size(Keys,1)
    if strncmp(Keys{i,1},'TTYPE',5)
        j=j+1;
        k=Keys{i,1};
        k=str2num(k(6:end));
        header{k}=Keys{i,2};
    end
end
if j~=a.BinaryTable(tableid).NFields
    error(['number of Keywords less than number of fields for file:',sprintf('\n'),...
        file,sprintf('\n'),num2str(j),',',num2str(a.BinaryTable(tableid).NFields)]);
end
fileinfo=a;