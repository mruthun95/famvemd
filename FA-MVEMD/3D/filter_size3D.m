
function Windows = filter_size3D(maxima_pos, minima_pos,type)
% Purpose:
% -To determine the window size for order statistics filtering of a signal.
% The determination of the window size is based on the work of Bhuiyan et
% al.
%
% Inputs: 
% -Two matrices of extrema positions
%
% Outputs:
% -Calculated value of the window size
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use delaunay triangulation to determine the nearest neighbours and hence
%the filter size

%processing d_max
max_nearest = zeros(length(maxima_pos),1);

try
TRI_max = delaunay(maxima_pos);
catch
    warning('Maxima points are collinear. Exiting without further iterations');
    Windows = [0, 0, 0, 0, 0, 0, 0];
    return
end
    

maxima_pos_x = maxima_pos(:,1);
maxima_pos_y = maxima_pos(:,2);
maxima_pos_z = maxima_pos(:,3);

%Calculating 6 edge lengths for each tetrahedron
e1 = sqrt( (maxima_pos_x(TRI_max(:,2))- maxima_pos_x(TRI_max(:,1))).^2 + (maxima_pos_y(TRI_max(:,2))- maxima_pos_y(TRI_max(:,1))).^2 + (maxima_pos_z(TRI_max(:,2))- maxima_pos_z(TRI_max(:,1))).^2 );
e2 = sqrt( (maxima_pos_x(TRI_max(:,3))- maxima_pos_x(TRI_max(:,1))).^2 + (maxima_pos_y(TRI_max(:,3))- maxima_pos_y(TRI_max(:,1))).^2 + (maxima_pos_z(TRI_max(:,3))- maxima_pos_z(TRI_max(:,1))).^2 );
e3 = sqrt( (maxima_pos_x(TRI_max(:,3))- maxima_pos_x(TRI_max(:,2))).^2 + (maxima_pos_y(TRI_max(:,3))- maxima_pos_y(TRI_max(:,2))).^2 + (maxima_pos_z(TRI_max(:,3))- maxima_pos_z(TRI_max(:,2))).^2 );
e4 = sqrt( (maxima_pos_x(TRI_max(:,4))- maxima_pos_x(TRI_max(:,1))).^2 + (maxima_pos_y(TRI_max(:,4))- maxima_pos_y(TRI_max(:,1))).^2 + (maxima_pos_z(TRI_max(:,4))- maxima_pos_z(TRI_max(:,1))).^2 );
e5 = sqrt( (maxima_pos_x(TRI_max(:,4))- maxima_pos_x(TRI_max(:,2))).^2 + (maxima_pos_y(TRI_max(:,4))- maxima_pos_y(TRI_max(:,2))).^2 + (maxima_pos_z(TRI_max(:,4))- maxima_pos_z(TRI_max(:,2))).^2 );
e6 = sqrt( (maxima_pos_x(TRI_max(:,4))- maxima_pos_x(TRI_max(:,3))).^2 + (maxima_pos_y(TRI_max(:,4))- maxima_pos_y(TRI_max(:,3))).^2 + (maxima_pos_z(TRI_max(:,4))- maxima_pos_z(TRI_max(:,3))).^2 );

%Calculating nearest neighbours for each maxima point
%Comparing tetrahedron edges associated with each vertex
em1 = min([e1, e2, e4],[],2); %Comparing edges 1, 2 and 4 (vertex 1)
em2 = min([e1, e3, e5],[],2); %Comparing edges 1, 3 and 5 (vertex 2)
em3 = min([e2, e3, e6],[],2); %Comparing edges 2, 3 and 6 (vertex 3)
em4 = min([e4, e5, e6],[],2); %Comparing edges 4, 5 and 6 (vertex 4)

e = [em1 ,em2, em3, em4];

%Making sure that the smallest edge associated with the each vertex is stored
%correctly
for i=1:length(em1)
    for j=1:4
        if max_nearest(TRI_max(i,j)) > e(i,j) || max_nearest(TRI_max(i,j)) == 0
            max_nearest(TRI_max(i,j)) = e(i,j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%processing d_min
min_nearest = zeros(length(minima_pos),1);

try
TRI_min = delaunay(minima_pos);
catch
    warning('Minima points are collinear. Exiting without further iterations');
    Windows = [0, 0, 0, 0, 0, 0, 0];
    return
end
minima_pos_x = minima_pos(:,1);
minima_pos_y = minima_pos(:,2);
minima_pos_z = minima_pos(:,3);

%Calculating 6 edge lengths for each tetrahedron
e1 = sqrt( (minima_pos_x(TRI_min(:,2))- minima_pos_x(TRI_min(:,1))).^2 + (minima_pos_y(TRI_min(:,2))- minima_pos_y(TRI_min(:,1))).^2 + (minima_pos_z(TRI_min(:,2))- minima_pos_z(TRI_min(:,1))).^2 );
e2 = sqrt( (minima_pos_x(TRI_min(:,3))- minima_pos_x(TRI_min(:,1))).^2 + (minima_pos_y(TRI_min(:,3))- minima_pos_y(TRI_min(:,1))).^2 + (minima_pos_z(TRI_min(:,3))- minima_pos_z(TRI_min(:,1))).^2 );
e3 = sqrt( (minima_pos_x(TRI_min(:,3))- minima_pos_x(TRI_min(:,2))).^2 + (minima_pos_y(TRI_min(:,3))- minima_pos_y(TRI_min(:,2))).^2 + (minima_pos_z(TRI_min(:,3))- minima_pos_z(TRI_min(:,2))).^2 );
e4 = sqrt( (minima_pos_x(TRI_min(:,4))- minima_pos_x(TRI_min(:,1))).^2 + (minima_pos_y(TRI_min(:,4))- minima_pos_y(TRI_min(:,1))).^2 + (minima_pos_z(TRI_min(:,4))- minima_pos_z(TRI_min(:,1))).^2 );
e5 = sqrt( (minima_pos_x(TRI_min(:,4))- minima_pos_x(TRI_min(:,2))).^2 + (minima_pos_y(TRI_min(:,4))- minima_pos_y(TRI_min(:,2))).^2 + (minima_pos_z(TRI_min(:,4))- minima_pos_z(TRI_min(:,2))).^2 );
e6 = sqrt( (minima_pos_x(TRI_min(:,4))- minima_pos_x(TRI_min(:,3))).^2 + (minima_pos_y(TRI_min(:,4))- minima_pos_y(TRI_min(:,3))).^2 + (minima_pos_z(TRI_min(:,4))- minima_pos_z(TRI_min(:,3))).^2 );

%Calculating nearest neighbours for each minima point
%Comparing tetrahedron edges associated with each vertex
emn1 = min([e1, e2, e4],[],2); %Comparing edges 1, 2 and 4 (vertex 1)
emn2 = min([e1, e3, e5],[],2); %Comparing edges 1, 3 and 5 (vertex 2)
emn3 = min([e2, e3, e6],[],2); %Comparing edges 2, 3 and 6 (vertex 3)
emn4 = min([e4, e5, e6],[],2); %Comparing edges 4, 5 and 6 (vertex 4)

e = [emn1 ,emn2, emn3, emn4];

%Making sure that the smallest edge associated with the each vertex is stored
%correctly
for i=1:length(emn1)
    for j=1:4
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