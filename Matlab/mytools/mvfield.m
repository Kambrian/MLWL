function s=mvfield(s,oldnames,newnames)
% rename fieldname from oldnames to newnames for structure s
%  oldnames and newnames can be cell arrays of the same size
%

if iscell(oldnames)
    for i=1:numel(oldnames)
        s.(newnames{i})=s.(oldnames{i});
    end
else
     s.(newnames)=s.(oldnames);
end
s=rmfield(s,oldnames);