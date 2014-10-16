function filename=gamaWL_fname(macros)
%construct the filename for gama

names=fieldnames(macros);
names=sort(names);
m=numel(names);
filename='gamawl_';
for i=1:m
    if ischar(macros.(names{i}))
    filename=[filename,macros.(names{i})];
    else
    filename=[filename,num2str(macros.(names{i}))];    
    end
end