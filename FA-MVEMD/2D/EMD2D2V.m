% Purpose: 
% -To perform EMD on 2 channels of 2 dimensional data
%
% Input: 
% - u: Signal 1
% - v: Signal 2
% - param
%   -nimfs: Number of IMFs to be extracted 
%   -tol: Sifting tolerance value
%   -type: type of window size to be used
%   -plot: 'on' to plot results, 'off' to hide IMF plots
%
% Output:
% - Results
%   - IMF (structure containing IMFs of all three signals)
%   - Residue (structure containing residue of all three signals)
%   - Windows (Window sizes (5 types) for each IMF)
%   - Sift_cnt (Number of sifting iterations for each signal)
%   - IO (Index of orthogonality for each signal)
%   - Error (Error of the decomposition for each signal)
%
% References:
%   [1] Bhuiyan et. al, 'Fast and Adaptive Bidimensional EmpiricalMode
%       Decomposition Using Order-Statistics Filter Based
%       Envelope Estimation',2008
%   
%   [2] FABEEMD (Matthew Koll, Dept. of Aerospace Engineering, University
%                of Illinois at Urbana-Champaign)
%
% 
% Written by Mruthun Thirumalaisamy
% Graduate Student
% Department of Aerospace Engineering
% University of Illinois at Urbana-Champaign
% May 11 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Results = EMD2D2V(u,v,param,varargin) 

%Reading signal characteristics
[Nx,Ny] = size(u); %Signal dimensions
B       = size(v); %Signal dimensions

%Preliminary checks
if ~isfield(param,'nimfs')
    param.nimfs = 10;
end

if ~isfield(param,'tol')
    param.tol = 0.05; % 0.1% of the minimum signal amplitude
end

if ~isfield(param,'type')
    param.type = 6;
end

if ~isfield(param,'plot')
    param.plot = 'off';
end

if(~all(ismember(param.type,[1,2,3,4,5,6,7])))
    error('Please enter a valid window size type')
end

if(~all([Nx,Ny]==B))
    error('Inconsistent dimensions between channels. Please check input data');
end
clearvars B

if(param.tol<=0.005)
   warning('Low sifting tolerance may cause oversifting');
   answer = questdlg('Would you like to continue?', ...
	'User set low sifting tolerance', ...
	'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            
        case 'No'
            return;
    end
end

%Initialisations
IMF.u = zeros(Nx, Ny ,param.nimfs); 
IMF.v = zeros(Nx, Ny ,param.nimfs);
Residue.u = u; Residue.v = v;

Windows = zeros(7,param.nimfs);

sift_cnt = zeros(1,param.nimfs);
imf = 1;
stopflag = 1;

    while(imf <= param.nimfs && stopflag)
        %Initialising intermediary IMFs
        H.u = Residue.u; H.v = Residue.v;

        sift_stop = 0; %flag to control sifting loop
        
        Combined = H.u/sqrt(2) + H.v/sqrt(2); %Combining two signals with equal weights  
        
        [maxima,minima] = Identify_max_min(Combined);  %Obtaining extrema of combined signal
        
        %Checking whether there are too few extrema in the IMF
        if (nnz(maxima) < 3 || nnz(minima) < 3)
            warning('Fewer than three extrema found in maxima map. Stopping now...');
            break;
        end
        
        %Window size determination by delaunay triangulation
        Windows(:,imf) = filter_size(maxima,minima,param.type);        
        w_sz = Windows(param.type,imf); %extracting window size chosen by input parameter
        
        if~(w_sz)
           warning('EMD2D3V has stopped because the Delaunay Triangulation could not be created (collinear points)'); 
           stopflag = 0; %#ok<NASGU>
           break;
        end
        
        %Begin sifting iteration
        while~(sift_stop)            
            sift_cnt(imf) = sift_cnt(imf) + 1; %Incrementing sift counter
            %Envelope Generation
            Env = OSF(H,w_sz);
            
            %padding
            Env = Pad_smooth(Env,w_sz);
           
            %Calculating mean envelope
            Env.u.med = (Env.u.maxs + Env.u.mins)./2;
            Env.v.med = (Env.v.maxs + Env.v.mins)./2;
            
            %Subtracting from residue
            H1.u = H.u - Env.u.med; H1.v = H.v - Env.v.med;       
            
            %Stop condition checks
            mse_u = immse(H1.u,H.u); mse_v = immse(H1.v,H.v);     
            if (mse_u<param.tol && mse_v<param.tol && sift_cnt(imf)~=1)
                sift_stop = 1;
            end
            
            H.u = H1.u; H.v = H1.v;            
        end
        
        %Storing IMFs
        IMF.u(:,:,imf) = H.u; IMF.v(:,:,imf) = H.v;

        %Subtracting from Residual Signals
        Residue.u = Residue.u - IMF.u(:,:,imf);
        Residue.v = Residue.v - IMF.v(:,:,imf);
        
        %Incrementing IMF counter
        imf = imf + 1;
        
    end
    
    %Checking for oversifting
    if(any(sift_cnt>=5*ones(size(sift_cnt))))
        warning('Decomposition may be oversifted. Checking if window size increases monotonically...');
        
        if( any (diff(Windows(param.type,:)) <= zeros(1,param.nimfs-1)) )
        warning('Filter window size does not increase monotonically')
        end
    end
    
    %Organising results
    Results.IMF = IMF;
    Results.Residue = Residue;
    Results.Windows = Windows;
    Results.Sifts = sift_cnt;
    
    %Error and orthogonality
    [Results.IO.u,Results.Error.u] = Orth_index(u,IMF.u,Residue.u);
    [Results.IO.v,Results.Error.v] = Orth_index(v,IMF.v,Residue.v);
    
    switch(param.plot)
        case 'on'
            Plot_results(u,v,Results,param)
    end
end

function [maxima,minima] = Identify_max_min(signal)
% Purpose:
% To identify the maxima and minima locations in a two or three dimensional array
mask = ones(3); mask(5) = 0; %Window size for the extrema detection fixed at 3x3 (Bhuiyan et.al)

B = ordfilt2(signal,8,mask);
C = ordfilt2(signal,1,mask);
maxima = signal >= B;
minima = signal <= C;
end

function Windows = filter_size(maxima_map, minima_map,type)
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

function Env = OSF(H,w_sz)
%Order statistics filtering to determine maximum and minmum envelopes
            Env.u.max = ordfilt2(H.u ,w_sz.^2, true(w_sz),'symmetric'); %Max envelope u
            Env.u.min = ordfilt2(H.u ,1, true(w_sz),'symmetric');       %Min envelope u
            
            Env.v.max = ordfilt2(H.v ,w_sz.^2, true(w_sz),'symmetric'); %Max envelope v
            Env.v.min = ordfilt2(H.v ,1, true(w_sz),'symmetric');       %Min envelope v            
end

function Env = Pad_smooth(Env,w_sz)
h = floor(w_sz/2);

%Padding
%u
Env.u.maxp = padarray(Env.u.max,[h h],'replicate');
Env.u.minp = padarray(Env.u.min,[h h],'replicate');
%v
Env.v.maxp = padarray(Env.v.max,[h h],'replicate');
Env.v.minp = padarray(Env.v.min,[h h],'replicate');

%Smoothing
%u
temp = movmean(Env.u.maxp,w_sz,2,'endpoints','discard');
Env.u.maxs = movmean(temp,w_sz,1,'endpoints','discard');
temp = movmean(Env.u.minp,w_sz,2,'endpoints','discard');
Env.u.mins = movmean(temp,w_sz,1,'endpoints','discard');
%v
temp = movmean(Env.v.maxp,w_sz,2,'endpoints','discard');
Env.v.maxs = movmean(temp,w_sz,1,'endpoints','discard');
temp = movmean(Env.v.minp,w_sz,2,'endpoints','discard');
Env.v.mins = movmean(temp,w_sz,1,'endpoints','discard');

end

function [IO,Error] = Orth_index(Signal,IMF,Residue)
% Purpose: 
% To calculate the index of orthogonality of a decomposition and its mean
% squared error

n_imf = size(IMF,3);
numerator = zeros(size(Signal));
I = sum(IMF,3) + Residue;

Error.map = (Signal-I)./Signal;
Error.global = immse(I,Signal);

for j = 1:n_imf
    for k = 1:n_imf
        if(j~=k)
           numerator = numerator + IMF(:,:,j).*IMF(:,:,k);
        end
    end
end

IO.map = numerator/sum(sum(Signal.^2));
IO.global = sum(sum(IO.map));
end

function Plot_results(u,v,Results,param)
% default plot attributes
set(groot,'defaultaxesfontname','times');
set(groot,'defaultaxesfontsize',12);
set(groot,'defaulttextInterpreter','latex');
set(groot,'defaultLineLineWidth',2);

Colour = redblue;

figure(1)   
        subplot(2,1,1)
        BIMF_plot(u,Colour,0,'Signal','u');
        subplot(2,1,2)
        BIMF_plot(v,Colour,0,'Signal','v');


    for i=1:param.nimfs
     figure(i+1)   
        subplot(2,1,1)
        BIMF_plot(Results.IMF.u(:,:,i),Colour,i,'IMF','u');
        subplot(2,1,2)
        BIMF_plot(Results.IMF.v(:,:,i),Colour,i,'IMF','v');
    end
    
    figure(i+2)
    subplot(2,1,1)
        BIMF_plot(Results.Residue.u,Colour,0,'Residue','u');
        subplot(2,1,2)
        BIMF_plot(Results.Residue.v,Colour,0,'Residue','v');
end

function BIMF_plot(signal,Colour,imf,name1,name2)
% %Masking wall data
% load('MASK_file','MASK');
% signal = MASK.*signal;

    imAlpha=ones(size(signal));
    imAlpha(isnan(signal))=0;    
    imagesc(signal,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    xlabel('$x$')
    ylabel('$y$')
    axis equal;
    axis tight;
    switch(name1)
        case 'IMF'
            title(sprintf('%s %d %s ',name1,imf,name2));
        case 'Residue'
            title(sprintf('%s %s ',name1,name2));
        case 'Signal'
            title(sprintf('%s %s ',name1,name2));
    end
    set(gca,'TickLabelInterpreter','latex')
    colormap(Colour);
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = 'u';
    set(colorTitleHandle ,'String',titleString,'Interpreter','latex','FontSize',14);
    set(hcb,'TickLabelInterpreter','latex');
end