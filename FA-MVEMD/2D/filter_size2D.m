
function Windows = filter_size2D(maxima_map, minima_map,type)
% Purpose:
% -To determine the window size for order statistics filtering of a signal.
% The determination of the window size is based on the work of Bhuiyan et
% al.
%
% Inputs: 
% -Two 2D extrema maps
%
% Outputs:
% -Calculated value of the window size
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use delaunay triangulation to determine the nearest neighbours and hence
%the filter size

%processing d_max
[maxima_pos_y,maxima_pos_x] = find(maxima_map);

max_nearest = zeros(length(maxima_pos_y),1);

try
TRI_max = delaunay([maxima_pos_x maxima_pos_y]);
catch
    warning('Maxima points are collinear. Exiting without further iterations');
    Windows = [0, 0, 0, 0, 0, 0, 0];
    return
end

%Calculating 3 edge lengths for each triangle
e1 = sqrt( (maxima_pos_x(TRI_max(:,2))- maxima_pos_x(TRI_max(:,1))).^2 + (maxima_pos_y(TRI_max(:,2))- maxima_pos_y(TRI_max(:,1))).^2 );
e2 = sqrt( (maxima_pos_x(TRI_max(:,3))- maxima_pos_x(TRI_max(:,1))).^2 + (maxima_pos_y(TRI_max(:,3))- maxima_pos_y(TRI_max(:,1))).^2 );
e3 = sqrt( (maxima_pos_x(TRI_max(:,3))- maxima_pos_x(TRI_max(:,2))).^2 + (maxima_pos_y(TRI_max(:,3))- maxima_pos_y(TRI_max(:,2))).^2 );

%Calculating nearest neighbours for each maxima point
%Comparing edges associated with each vertex
em1 = min([e2, e1],[],2); %Comparing edges 2 and 1 (vertex 1)
em2 = min([e3, e1],[],2); %Comparing edges 3 and 1 (vertex 2)
em3 = min([e3, e2],[],2); %Comparing edges 3 and 2 (vertex 3)

e = [em1 ,em2, em3]; %Smaller edge for each vertex in each triangle (since one vertex is associated with two edges in a triangle)

%Making sure that the smallest edge associated with the each vertex is stored
%correctly
for i=1:length(em1)
    for j=1:3
        if max_nearest(TRI_max(i,j)) > e(i,j) || max_nearest(TRI_max(i,j)) == 0
            max_nearest(TRI_max(i,j)) = e(i,j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%processing d_min

[minima_pos_y,minima_pos_x] = find(minima_map);
min_nearest = zeros(length(minima_pos_y),1);

try
TRI_min = delaunay([minima_pos_x minima_pos_y]);
catch
    warning('Minima points are collinear. Exiting without further iterations');
    Windows = [0, 0, 0, 0, 0, 0, 0];
    return
end

%Calculating 3 neighbour distances for each minima point
e1 = sqrt( (minima_pos_x(TRI_min(:,2))- minima_pos_x(TRI_min(:,1))).^2 + (minima_pos_y(TRI_min(:,2))- minima_pos_y(TRI_min(:,1))).^2 );
e2 = sqrt( (minima_pos_x(TRI_min(:,3))- minima_pos_x(TRI_min(:,1))).^2 + (minima_pos_y(TRI_min(:,3))- minima_pos_y(TRI_min(:,1))).^2 );
e3 = sqrt( (minima_pos_x(TRI_min(:,3))- minima_pos_x(TRI_min(:,2))).^2 + (minima_pos_y(TRI_min(:,3))- minima_pos_y(TRI_min(:,2))).^2 );

%Calculating nearest neighbours for each maxima point

%Comparing triangle edges associated with each vertex
emn1 = min([e2, e1],[],2); %Comparing edges 2 and 1 (vertex 1)
emn2 = min([e3, e1],[],2); %Comparing edges 3 and 1 (vertex 2)
emn3 = min([e3, e2],[],2); %Comparing edges 3 and 2 (vertex 3)

e = [emn1 ,emn2, emn3]; %Smaller edge for each vertex in each triangle (since one vertex is associated with two edges in a triangle)

%Making sure that the smallest edge associated with the each vertex is stored
%correctly
for i=1:length(emn1)
    for j=1:3
        if min_nearest(TRI_min(i,j)) > e(i,j) || min_nearest(TRI_min(i,j)) == 0
            min_nearest(TRI_min(i,j)) = e(i,j);
        end
    end
end

%Window size calculations

d1 = min( min(max_nearest) , min(min_nearest) );
d2 = max( min(max_nearest) , min(min_nearest) );
d3 = min( max(max_nearest) , max(min_nearest) );
d4 = max( max(max_nearest) , max(min_nearest) );
d5 = (d1+d2+d3+d4)/4 ;
d6 = median([min_nearest; max_nearest]);
d7 = mode([min_nearest; max_nearest]);

Windows = [d1, d2, d3, d4, d5, d6, d7];

%making sure w_size is an odd integer
Windows = 2*(floor(Windows./2))+1;
         
if(Windows(type)<3)
    warning('WARNING: Calculated Window size less than 3');
    warning('Overriding calculated value and setting window size = 3');
    Windows(type) = 3;
end
end
