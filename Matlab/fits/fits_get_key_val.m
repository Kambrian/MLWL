function val=fits_get_key_val( names,Keywords )
% val=fits_get_key_val( names,Key )
% to get the value of keys from fits keywords data for names
% input: names, cell array of names to query
%           Keywords, cell array of name-val-comments data
% output: val, values corresponding to names

if ~iscell(names)
        warning('input is not cell-array, assuming each row to be one keyname');
        n=size(names,1);
        names=mat2cell(names,ones(1,n),size(names,2));
else
n=numel(names);
end
val=cell(size(names));

m=size(Keywords,1);
for i=1:n
    for j=1:m
        if strcmp(names{i},Keywords{j,1})
            val{i}=Keywords{j,2};
            break;
        end
    end
end

end

