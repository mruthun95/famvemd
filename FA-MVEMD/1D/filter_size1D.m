function window_size = filter_size1D(imax, imin)
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

edge_len_max = diff(imax);
edge_len_min = diff(imin);
    
    
        %Window size calculations
        
        d1 = min( min(edge_len_max) , min(edge_len_min) );
        d2 = max( min(edge_len_max) , min(edge_len_min) );
        d3 = min( max(edge_len_max) , max(edge_len_min) );
        d4 = max( max(edge_len_max) , max(edge_len_min) );
        d5 = (d1+d2+d3+d4)/4 ;
 
        %making sure w_size is an odd integer
         window_size = 2*(floor(d4/2))+1;
        
if(window_size<3)
    warning('WARNING: Calculated Window size less than 3\n ');
    warning('Overriding calculated value and setting window size = 3');
    window_size = 3;
end
         
    

        
        