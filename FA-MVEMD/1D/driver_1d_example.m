%cleanup
clearvars
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
param.n_grid_1    = 500;  %Grid size in dimension 1
param.nimfs       = 2;    %Maximum number of IMFs that can be stored
param.type        = 5; %type of window size
param.tol         = 0.05; %sifting tolerance
param.plot        = 'on'; %plots on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = linspace(0,6*pi,param.n_grid_1);

u = 2.5*cos(t); %Simple Signal

v = 2.5*cos(5*t); %Simple Signal

w = u+v;

Results = EMD1D3V(u',v',w',t,param); %Make sure signals are column vectors while being passed to EMD function



