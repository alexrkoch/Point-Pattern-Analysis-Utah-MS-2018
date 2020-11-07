% VerticeInterpolator.m
% Developed by Alex Koch, contact at alexrkoch@gmail.com
% Repository at https://github.com/alexrkoch

% Interpolate vertices of a set of polygons so there is one for every 
% meter or decimeter of X values. This is necessary in order to use
% the width, thickness, NTG, and amalgamation scripts of this thesis.


clear;

%Load polygon-bearing shapefiles or structures.
shp = shaperead('CC5_Flattened_2D.shp'); % input original shapefile or
    % .mat structure.

saveName = ('CC5_Shapefile_mInterp.mat'); % choose name to save to.

%Set Interpolation Interval
interp = 1; % for every meter use '1'. for every decimeter use '0.1'
interpFactor = 1; % for meter use '1'. for decimeter use '100'

%% Loop through each polygon and perform the workflow

% Remove the NaN from last cell of each array. When importing 
% shapefiles from ArcGIS, the last cell of each polygon's X-Y vertice
% dataset is NaN. When this is the case, MATLAB won't calculate area,
% which we need for NTG. So just check your file structure first, if 
% there is no NaN at the end, then comment this following loop out.
for a = 1:length(shp)
    lastX = length(shp(a).X);
    shp(a).X(lastX) = [];
    shp(a).Y(lastX) = [];
end

CC5shp = shp; %Change to current variable name

for a = 1:length(shp) % search each polygon
newXvertices = []; %will be replaced on the next polygon
newYvertices = [];
for b = 1:length(shp(a).X)
        newXvertices = horzcat(newXvertices,shp(a).X(b)); 
            % add X-val b vertice to new structure.            
        newYvertices = horzcat(newYvertices,shp(a).Y(b));  
            % add Y-val bvertice to new structure.
        if b <= length(shp(a).X)-2 % So the loop stops when it gets to 
            % the penultimate vertice (last vertice is a copy of 
            % the first vertice)
        if (shp(a).X(b) < shp(a).X(b+1)) %if vertice X values are
            %increasing with increasing index
            c = b+1 ;
            polyStart = ceil((shp(a).X(b))*(interp*interpFactor))/ ...
                (interp*interpFactor); %find the first even meter or  
                % decimeter greater than vertice b.
            polyEnd = floor((shp(a).X(c))*(interp*interpFactor))/ ... 
                (interp*interpFactor); %find the first even meter or 
                % decimeter less than vertice c.
            m = ((shp(a).Y(b))-(shp(a).Y(c)))/((shp(a).X(b))- ...
                (shp(a).X(c))); % slope = (Y1-Y2) / (X1-X2)
            f = (polyEnd-polyStart)/interp ;
                for d = 0:f % d are the new regularly spaced X values.
                    g = polyStart + (d*interp); % new X value
                    e = m*(g-(shp(a).X(b)))+(shp(a).Y(b));
                        % calculate new Y value
                    newXvertices = horzcat(newXvertices,g);
                        % add new x and y (below) to new structure
                    newYvertices = horzcat(newYvertices,e);
                end
        elseif (shp(a).X(b) > shp(a).X(b+1))  %if vertice X values are 
            % decreasing with increasing index
            c = b+1 ;                
            polyEnd = floor((shp(a).X(b))*(interp*interpFactor))/ ... 
                (interp*interpFactor); %find the first even meter or
                % decimeter greater than vertice b.
            polyStart = ceil((shp(a).X(c))*(interp*interpFactor))/ ...
                (interp*interpFactor); %find the first even meter or
                % decimeter less than vertice c.
            m = ((shp(a).Y(b))-(shp(a).Y(c)))/((shp(a).X(b))- ...
                (shp(a).X(c))); % slope = (Y1-Y2) / (X1-X2)
            f = (polyEnd-polyStart)/interp ;
                for d = 0:f % d are the new regularly  
                    % spaced X values.
                    g = polyStart + (d*interp); % new X value
                    e = m*(g-(shp(a).X(b)))+(shp(a).Y(b));
                        % calculate new Y value
                    newXvertices = horzcat(newXvertices,g);
                        % add new x and y (below) to new structure
                    newYvertices = horzcat(newYvertices,e);
                end
        end % end if loop for vertices b and c.
        end % end index control
    newShp(a).X = newXvertices;
    newShp(a).Y = newYvertices;

end % end b vertice for
end % end for 'a' shape

% Add the new Vertices to the final line file.
for a = 1:length(shp)
    CC5shp(a).X = newShp(a).X;   %Change to current variable name
    CC5shp(a).Y = newShp(a).Y;   %Change to current variable name
end

save(saveName,'CC5shp');
