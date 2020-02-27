function Windows = filter_size1D(imax, imin, type)
%
% Purpose:
% -To determine the window size for order statistics filtering of a signal.
% The determination of the window size is based on the work of Bhuiyan et
% al.
%
% Inputs:
% -Two 1D extrema maps
%
% Outputs:
% -Calculated value of the window size
%
% Written by Mruthun Thirumalaisamy
% Graduate Student
% Department of Aerospace Engineering
% University of Illinois at Urbana-Champaign
% March 30 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edge_len_max = diff(sort(imax));
edge_len_min = diff(sort(imin));


%Window size calculations

d1 = min( min(edge_len_max) , min(edge_len_min) );
d2 = max( min(edge_len_max) , min(edge_len_min) );
d3 = min( max(edge_len_max) , max(edge_len_min) );
d4 = max( max(edge_len_max) , max(edge_len_min) );
d5 = (d1+d2+d3+d4)/4 ;
d6 = median([edge_len_max; edge_len_min]);
d7 = mode([edge_len_max; edge_len_min]);

Windows = [d1, d2, d3, d4, d5, d6, d7];

%making sure w_size is an odd integer
Windows = 2*(floor(Windows./2))+1;

if(Windows(type)<3)
    warning('WARNING: Calculated Window size less than 3');
    warning('Overriding calculated value and setting window size = 3');
    Windows(type) = 3;
end

end