function table=fits_load_bintable(file,cols,tableid, format) %#ok<STOUT>
%function table=fits_load_bintable(file,cols,tableid,format)
%
% to load raw binary table contents from fits file
%
% input: file: fits file to read
%        cols: optional, a vector listing the indexes of fields to load;
%                   default=0: load all cols.
%        tableid: optional, index of the bintable to load in case there're
%                    more than one bintalbes available; default=1.
%       format: optional, boolen, 0: read the data in raw
%                    format; 1: the data would be scaled. default=0
% output: table: optional, struct containing all the fields read;
%            if output is not explicitly specified, load the fields directly
%            into caller's memory space.
% example: table=fits_load_bintable('sdss_01.fits',[1,3,5]) will load the
% 1st,3rd and fifth fields from the bintable as subfields into table
% without assign output, fits_load_bintable('sdss_01.fits',[1,3,5]) will
% create the [1,3,5] fields as individual variables
% 
% ==Jiaxin Han, 11/02/2010, Durham ==
% Revision: 06/16/2011, Durham


if nargin<3, tableid=1; end
if nargin<4, format=0; end

[header,a]=fits_bintable_keys(file,tableid);
if sum(a.BinaryTable(tableid).Intercept~=0)||sum(a.BinaryTable(tableid).Slope~=1)
    warn('Binary table has scale, need to be scale manually');
end

if format
    data=fitsread(file,'bintable',tableid,'info',a);
else
    data=fitsread(file,'bintable',tableid,'raw','info',a);
end

if nargin<2||cols(1)==0
    cols=1:numel(header);
end
%create variables
for i=1:numel(cols)
    j=cols(i);
    if ~isvarname(header{j})
        header{j}=genvarname(header{j});
    end
    if nargout
%     eval(['table.',header{j},'=data{j};']);
    table.(header{j})=data{j};
    else
    assignin('caller',header{j},data{j});
    end
end
    
