%cleanup
clearvars
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
param.n_grid_1    = 30;  %Grid size in dimension 1
param.n_grid_2    = 30;  %Grid size in dimension 2
param.n_grid_3    = 30;  %Grid size in dimension 3
param.nimfs       = 3;    %Maximum number of IMFs that can be stored
param.type        = 5;    %type of window size
param.tol         = 0.05; %sifting tolerance
param.nslice      = 10;   %number of slices in volume plot
param.plot        = 'on'; %toggle plotting

%USE FOR EMD3D3V_parallel_var
% param.xend        = 5; %x-domain length in physical units 
% param.yend        = 5; %y-domain length in physical units
% param.zend        = 5; %z-domain length in physical units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = linspace(0,30,param.n_grid_1);
Y = linspace(0,15,param.n_grid_2);
Z = linspace(0,15,param.n_grid_3);

[x,y,z] = meshgrid(Y,X,Z); 

uc1 = (sin(3.5*x)+sin(3.5*y)+sin(3.5*z));  uc2 = (sin(x)+sin(y)+sin(z)); ru  = (cos(x/24)+cos(y/24)+cos(z/24));
u   =  uc1 + uc2 + ru; %Signal 'u'

vc1 = (sin(x)+sin(y)+sin(z));  vc2 = (sin(3.5*x)+sin(3.5*y)+sin(3.5*z)); rv  = (cos(x/24)+cos(y/24)+cos(z/24));
v   =  vc1 ; %Signal 'v'

wc1 = (cos(x/24)+cos(y/24)+cos(z/24));  wc2 = (sin(x/5)+sin(y/5)+sin(z/5)); rw = (sin(x)+sin(y)+sin(z)) ;
w   =  wc1 ;%+ wc2 + rw; %Signal 'w'

Results = EMD3D3V(u,v,w,param);



