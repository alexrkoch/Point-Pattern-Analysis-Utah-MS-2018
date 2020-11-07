% This script is the same as from Cell.m from Benhallam (2015), except
% with an added field 'z' which is used to store a property for a cell.

% This function is called on in the CellColorMaps.m script.

classdef colorCell
    properties
        x;
        y;
        width;
        height;
        z;
    end
    
    methods
        %constructor
        function obj = colorCell(x,y,width,height,z)
            if nargin ~= 0
                obj.x = x;
                obj.y = y;
                obj.width = width;
                obj.height = height;
                obj.z = z;
            end
        end
    end
end
