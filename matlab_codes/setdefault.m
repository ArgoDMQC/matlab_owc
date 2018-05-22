function fl = setdefault(fl,field,value)
% set a default value in the field of a structure.
%fl is a struct array
%field is a string variable containg the name of a field.

if isfield(fl,field)
    return
end
eval(['fl.' field ' = value;']);