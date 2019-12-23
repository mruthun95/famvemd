%cleanup
clearvars
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
param.n_grid_1    = 256;  %Grid size in dimension 1
param.n_grid_2    = 256;  %Grid size in dimension 2
param.nimfs       = 3;    %Maximum number of IMFs that can be stored
param.type        = 6;    %type of window size
param.tol         = 0.01; %sifting tolerance
param.plot        = 'on'; %toggle plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = linspace(0,30,param.n_grid_1);
Y = linspace(0,30,param.n_grid_2);

[x,y] = meshgrid(Y,X); 

uc1 = (sin(4*x)+sin(4*y));  uc2 = (sin(x)+sin(y)); ru  = (cos(x/24)+cos(y/24));
u   =  uc1 + uc2 + ru; %Signal 'u'

vc1 = (sin(x)+sin(y));  vc2 = (sin(3.5*x)+sin(3.5*y)); rv  = (cos(x/24)+cos(y/24));
v   =  vc1 ; %Signal 'v'

wc1 = (cos(x/24)+cos(y/24));  wc2 = (sin(x/5)+sin(y/5)); rw = (sin(x)+sin(y)) ;
w   =  wc2 ;%+ wc2 + rw; %Signal 'w'

Results = EMD2D3V(u,v,w,param);



