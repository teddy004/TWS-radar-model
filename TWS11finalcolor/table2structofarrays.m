function [ outStruct ] = table2structofarrays( inTable )
%TABLE2STRUCTOFARRAYS Convert a table to a struct of arrays.
%   Usage:  outStruct = TABLE2STRUCTOFARRAYS( inTable )
%
% Convert a table with M rows and N variables to a struct with N fields, 
% each of which contains a 1-dimensional array of length M and where the
% field names in the struct are the same as the variable names in the
% table.
%
% NOTE: There ia a built-in function TABLE2STRUCT which converts a table to
% an array of structs.  However, there are HUGE performance advantages of
% a struct of arrays over an array of structs.

% Make sure the input really is a table
if ~isa(inTable, 'table')
    error('Error.  Input to function %s must be a table, not a %s', mfilename, class(inTable))
end

% Create an empty struct with no fields
outStruct = struct;

% Iterate through all of the variables in the table
for varNum=1:width(inTable)
    % Get the variable name as a cell array with 1 element
    varNameCell = inTable.Properties.VariableNames(varNum);
    
    % Extract the variable name as a string
    varName = varNameCell{1};
    
    % Add a new field to the struct containing the data for this variable
    outStruct = setfield(outStruct, varName, inTable.(varNum));
end

% If the table has explicitly defined row names, then add a field for these
% NOTE: This is done last to preserve variable number order
if ~isempty(inTable.Properties.RowNames)
    outStruct = setfield(outStruct, 'RowNames', inTable.Properties.RowNames)
end

end

