function [maxima,minima] = Identify_max_min(signal)
% Purpose:
% To identify the maxima and minima locations in a two or three dimensional array
mask = ones(3); mask(5) = 0; %Window size for the extrema detection fixed at 3x3 (Bhuiyan et.al)

B = ordfilt2(signal,8,mask);
C = ordfilt2(signal,1,mask);
maxima = signal >= B;
minima = signal <= C;
end