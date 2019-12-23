% Purpose: 
% -To perform EMD on 3 channels of 3 dimensional data
%
% Input: 
% - u: Signal 1
% - v: Signal 2
% - w: Signal 3
% - param.
%   -nimfs: Number of IMFs to be extracted 
%   -tol: Sifting tolerance value
%   -type: type of window size to be used
%   -plot: 'on' to plot results, default hides IMF plots
%   -nslice: number of slices in volume plot
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
%
% 
% Written by Mruthun Thirumalaisamy
% Graduate Student
% Department of Aerospace Engineering
% University of Illinois at Urbana-Champaign
% May 16 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Results = EMD3D3V(u,v,w,param) 

%Reading signal characteristics
[Nx,Ny,Nz] = size(u); %Signal dimensions
B          = size(v); %Signal dimensions
C          = size(w); %Signal dimensions

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

if(~all([Nx,Ny,Nz]==B) ||  ~all([Nx,Ny,Nz]==C))
    error('Inconsistent dimensions between channels. Please check input data');
end
clearvars B C

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
IMF.u = zeros(Nx, Ny, Nz, param.nimfs); 
IMF.v = zeros(Nx, Ny, Nz, param.nimfs);
IMF.w = zeros(Nx, Ny, Nz, param.nimfs);
Residue.u = u; Residue.v = v; Residue.w = w;

Windows = zeros(7,param.nimfs);

sift_cnt = zeros(1,param.nimfs);
imf = 1;
stopflag = 1;

    while(imf <= param.nimfs && stopflag)
        %Initialising intermediary IMFs
        H.u = Residue.u; H.v = Residue.v; H.w = Residue.w;

        sift_stop = 0; %flag to control sifting loop
        
        Combined = H.u/sqrt(3) + H.v/sqrt(3) + H.w/sqrt(3); %Combining three signals with equal weights        
        [Maxima,MaxPos,Minima,MinPos] = MinimaMaxima3D(Combined,1,1,[],[]);  %Obtaining extrema of combined signal
        
        %Checking whether there are too few extrema in the IMF
        if (nnz(Maxima) < 3 || nnz(Minima) < 3)
            warning('Fewer than three extrema found in extrema map. Stopping now...');
            break;
        end
        
        %Window size determination by delaunay triangulation
        Windows(:,imf) = filter_size(MaxPos,MinPos,param.type);        
        w_sz = Windows(param.type,imf); %extracting window size chosen by input parameter
        
        if~(w_sz)
           warning('EMD3D3V has stopped because the Delaunay Triangulation could not be created (collinear points)'); 
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
            Env.w.med = (Env.w.maxs + Env.w.mins)./2;
            
            %Subtracting from residue
            H1.u = H.u - Env.u.med; H1.v = H.v - Env.v.med; H1.w = H.w - Env.w.med;         
             
            %Stop condition checks
            mse_u = immse(H1.u,H.u); mse_v = immse(H1.v,H.v); mse_w = immse(H1.w,H.w);         
            if (mse_u<param.tol && mse_v<param.tol && mse_w<param.tol && sift_cnt(imf)~=1)
                sift_stop = 1;
            end
            
            H.u = H1.u; H.v = H1.v; H.w = H1.w;                
        end
        
        %Storing IMFs
        IMF.u(:,:,:,imf) = H.u; IMF.v(:,:,:,imf) = H.v; IMF.w(:,:,:,imf) = H.w;

        %Subtracting from Residual Signals
        Residue.u = Residue.u - IMF.u(:,:,:,imf);
        Residue.v = Residue.v - IMF.v(:,:,:,imf);
        Residue.w = Residue.w - IMF.w(:,:,:,imf);  
        
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
    Results.windowtype = param.type;
    Results.Residue = Residue;
    Results.Windows = Windows;
    Results.Sifts = sift_cnt;
    
    %Error and orthogonality
    [Results.IO.u,Results.Error.u] = Orth_index(u,IMF.u,Residue.u);
    [Results.IO.v,Results.Error.v] = Orth_index(v,IMF.v,Residue.v);
    [Results.IO.w,Results.Error.w] = Orth_index(w,IMF.w,Residue.w);
    
    switch(param.plot)
        case 'on'
            Plot_results(u,v,w,Results,param)
    end
end

function Windows = filter_size(maxima_pos, minima_pos,type)
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

function Env = OSF(H,w_sz)
%Order statistics filtering to determine maximum and minmum envelopes
            Env.u.max = Separable_ordfilt3(H.u, 'max', w_sz); %Max envelope u
            Env.u.min = Separable_ordfilt3(H.u, 'min', w_sz); %Min envelope u
             
            Env.v.max = Separable_ordfilt3(H.v, 'max', w_sz); %Max envelope v
            Env.v.min = Separable_ordfilt3(H.v, 'min', w_sz); %Min envelope v
            
            Env.w.max = Separable_ordfilt3(H.w, 'max', w_sz); %Max envelope w
            Env.w.min = Separable_ordfilt3(H.w, 'min', w_sz); %Min envelope w
            
            function Signal = Separable_ordfilt3(Signal, order, w_sz)
                % Purpose:
                % -To perform separable order statistics filtering of 3D
                % signals 
                % -Boundary condition is always symmetric
               
                [X,Y,Z] = size(Signal);
                
                %Separable Filtering
                %First Dimension (X)
                for k = 1:Z
                    for j = 1:Y
                        Signal(:,j,k) = Ordfilt1(Signal(:,j,k),order,w_sz);
                    end
                end
                
                %Second Dimension (Y)
                for k = 1:Z
                    for i = 1:X
                        Signal(i,:,k) = Ordfilt1(Signal(i,:,k),order,w_sz);
                    end
                end
                
                %Third Dimension (Z)
                for j = 1:Y
                    for i = 1:X
                        Signal(i,j,:) = Ordfilt1(Signal(i,j,:),order,w_sz);
                    end
                end
                
                function f_signal = Ordfilt1(signal,order,window_size)
                    
                    %1-D Rank order filter function
                    
                    %Pre-processing
                    [a,b,c] = size(signal);           %Original signal size
                    signal  = squeeze(signal);        %Removing the singleton dimensions
                    L       = length(signal);         %Length of the signal
                    signal  = reshape(signal, [L,1]); %Ensure that the processed signal is always a column vector
                    
                    r = (window_size-1)/2;
                    
                    %Padding boundaries
                    x = [flip(signal(1:r)); signal ;flip(signal(end-(r-1):end))];
                    
                    [M,~] = size(x);
                    y = zeros(size(x));
                                            
                    switch order
                        case 'max'
                            for m = 1+r:M-r
                                % Extract a window of size (2r+1) around (m)
                                temp = x((m-r):(m+r));
                                w = sort(temp);
                                y(m) = w(end); % Select the greatest element
                            end
                        case 'min'
                            for m = 1+r:M-r
                                % Extract a window of size (2r+1) around (m)
                                temp = x((m-r):(m+r));
                                w = sort(temp);
                                y(m) = w(1); % Select the smallest element
                            end
                        otherwise
                            error('No such filering operation defined')
                    end
                    
                    f_signal = y(1+r:end-r);
                    
                    f_signal = reshape(f_signal,[a,b,c]); %Restoring Signal size   
                end          
            end
end

function Env = Pad_smooth(Env,w_sz)
h = floor(w_sz/2);

%Padding
%u
temp = padarray(Env.u.max,[h h],'replicate');
temp1 = permute(temp,[3 2 1]); %interchanging dimensions
temp = padarray(temp1,[h 0],'replicate');
Env.u.maxp = permute(temp,[3 2 1]); %restoring dimensions

temp = padarray(Env.u.min,[h h],'replicate');
temp1 = permute(temp,[3 2 1]); %interchanging dimensions
temp = padarray(temp1,[h 0],'replicate');
Env.u.minp = permute(temp,[3 2 1]); %restoring dimensions
%v
temp = padarray(Env.v.max,[h h],'replicate');
temp1 = permute(temp,[3 2 1]); %interchanging dimensions
temp = padarray(temp1,[h 0],'replicate');
Env.v.maxp = permute(temp,[3 2 1]); %restoring dimensions

temp = padarray(Env.v.min,[h h],'replicate');
temp1 = permute(temp,[3 2 1]); %interchanging dimensions
temp = padarray(temp1,[h 0],'replicate');
Env.v.minp = permute(temp,[3 2 1]); %restoring dimensions
%w
temp = padarray(Env.w.max,[h h],'replicate');
temp1 = permute(temp,[3 2 1]); %interchanging dimensions
temp = padarray(temp1,[h 0],'replicate');
Env.w.maxp = permute(temp,[3 2 1]); %restoring dimensions

temp = padarray(Env.w.min,[h h],'replicate');
temp1 = permute(temp,[3 2 1]); %interchanging dimensions
temp = padarray(temp1,[h 0],'replicate');
Env.w.minp = permute(temp,[3 2 1]); %restoring dimensions

%Smoothing
%u
temp1 = movmean(Env.u.maxp,w_sz,3,'endpoints','discard');
temp2 = movmean(temp1,w_sz,2,'endpoints','discard');
Env.u.maxs = movmean(temp2,w_sz,1,'endpoints','discard');
temp1 = movmean(Env.u.minp,w_sz,3,'endpoints','discard');
temp2 = movmean(temp1,w_sz,2,'endpoints','discard');
Env.u.mins = movmean(temp2,w_sz,1,'endpoints','discard');
%v
temp1 = movmean(Env.v.maxp,w_sz,3,'endpoints','discard');
temp2 = movmean(temp1,w_sz,2,'endpoints','discard');
Env.v.maxs = movmean(temp2,w_sz,1,'endpoints','discard');
temp1 = movmean(Env.v.minp,w_sz,3,'endpoints','discard');
temp2 = movmean(temp1,w_sz,2,'endpoints','discard');
Env.v.mins = movmean(temp2,w_sz,1,'endpoints','discard');
%w
temp1 = movmean(Env.w.maxp,w_sz,3,'endpoints','discard');
temp2 = movmean(temp1,w_sz,2,'endpoints','discard');
Env.w.maxs = movmean(temp2,w_sz,1,'endpoints','discard');
temp1 = movmean(Env.w.minp,w_sz,3,'endpoints','discard');
temp2 = movmean(temp1,w_sz,2,'endpoints','discard');
Env.w.mins = movmean(temp2,w_sz,1,'endpoints','discard');

end

function [IO,Error] = Orth_index(Signal,IMF,Residue)
% Purpose: 
% To calculate the index of orthogonality of a decomposition and its mean
% squared error

n_imf = size(IMF,4);
numerator = zeros(size(Signal));
I = sum(IMF,4) + Residue;

Error.map = (Signal-I)./Signal;
Error.global = immse(I,Signal);

for j = 1:n_imf
    for k = 1:n_imf
        if(j~=k)
           numerator = numerator + IMF(:,:,:,j).*IMF(:,:,:,k);
        end
    end
end

IO.map = numerator/sum(sum(sum(Signal.^2))); %wrong
IO.global = sum(sum(sum(IO.map)));
end

function Plot_results(u,v,w,Results,param)
% default plot attributes
set(groot,'defaultaxesfontname','times');
set(groot,'defaultaxesfontsize',12);
set(groot,'defaulttextInterpreter','latex');
set(groot,'defaultLineLineWidth',2);

Colour = redblue();
nslice  = param.nslice;

figure(1)   
        subplot(1,3,1)
        TIMF_plot(u,Colour,nslice,0,'Signal','u');
        subplot(1,3,2)
        TIMF_plot(v,Colour,nslice,0,'Signal','v');
        subplot(1,3,3)
        TIMF_plot(w,Colour,nslice,0,'Signal','w');


    for i=1:param.nimfs
     figure(i+1)   
        subplot(1,3,1)
        TIMF_plot(Results.IMF.u(:,:,:,i),Colour,nslice,i,'IMF','u');
        subplot(1,3,2)
        TIMF_plot(Results.IMF.v(:,:,:,i),Colour,nslice,i,'IMF','v');
        subplot(1,3,3)
        TIMF_plot(Results.IMF.w(:,:,:,i),Colour,nslice,i,'IMF','w');
    end
    
    figure(i+2)
    subplot(1,3,1)
        TIMF_plot(Results.Residue.u,Colour,nslice,0,'Residue','u');
        subplot(1,3,2)
        TIMF_plot(Results.Residue.v,Colour,nslice,0,'Residue','v');
        subplot(1,3,3)
        TIMF_plot(Results.Residue.w,Colour,nslice,0,'Residue','w');
end

function TIMF_plot(signal,Colour,nslice,imf,name1,name2)    

    [Nx, Ny, Nz] = size(signal);

    xslice = linspace(1,Nx,nslice);
    yslice = linspace(1,Ny,nslice);
    zslice = linspace(1,Nz,nslice);
    volume = slice(signal,xslice,yslice,zslice);
    axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    set(gca,'TickLabelInterpreter','latex')
    switch(name1)
        case 'IMF'
            title(sprintf('%s %d %s',name1,imf,name2));
        case 'Signal'
            title(sprintf('%s %s',name1,name2));
        case 'Residue'
            title(sprintf('%s %s',name1,name2));
    end
    colorbar;
    set(volume,'EdgeColor','none',...
        'FaceColor','interp',...
        'FaceAlpha','interp')
    alpha('color')
    view(30,30);
    alphamap('rampup')
    alphamap('decrease',.1)
    colormap(Colour);
    caxis([-3 3]);
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = sprintf('%s',name2);
    set(colorTitleHandle ,'String',titleString,'Interpreter','latex','FontSize',14);
    set(hcb,'TickLabelInterpreter','latex');
    
end