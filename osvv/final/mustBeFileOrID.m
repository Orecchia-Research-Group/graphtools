function [] = mustBeFileOrID(filename, read, write)
%MUSTBEFILE Checks that argument is either a valid file or 
% stdin, stdout, stderr (0, 1, 2) respectively
arguments
    filename
    read (1, 1) {mustBeNumericOrLogical} = true;
    write (1, 1) {mustBeNumericOrLogical} = false;
end

if isa(filename, 'char')
    mustBeFile(filename);
    [~, attribs] = fileattrib(filename);
    if read && ~attribs.UserRead
        msgType = sprintf('File %s needs to be readable, but is not', filename);
        eidType = 'mustBeFile:notReadable';
        error(eidType, msgType);
    end
    if write && ~attribs.UserWrite
        msgType = sprintf('File %s needs to be writable, but is not', filename);
        eidType = 'mustBeFile:notWriteable';
        error(eidType, msgType);
    end
else
    mustBeNumeric(filename)
    if read && filename ~= 0
        msgType = 'Argument file stream needs to be stdin (0) to be readable.';
        eidType = 'mustBeFile:notReadable';
        error(eidType, msgType);
    elseif write && (filename ~= 1 && filename ~= 2)
        msgType = 'Argument file stream needs to be stdout (1) or stderr (2) to be readable.';
        eidType = 'mustBeFile:notWriteable';
        error(eidType, msgType);
    end
end


