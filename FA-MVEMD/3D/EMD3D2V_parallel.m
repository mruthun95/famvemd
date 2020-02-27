% Purpose: 
% -To perform EMD on 2 channels of 3 dimensional data
%
% Input: 
% - u: Signal 1
% - v: Signal 2
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
% May 16 2018 (Modified: Dec 14 2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Results = EMD3D2V_parallel(u,v,param) 

%Reading signal characteristics
[Nx,Ny,Nz] = size(u); %Signal dimensions
B          = size(v); %Signal dimensions

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

if(~all([Nx,Ny,Nz]==B))
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

if ~isfield(param,'plot')
    param.plot = 'off';
end

%Initialisations
IMF.u = zeros(Nx, Ny, Nz, param.nimfs); 
IMF.v = zeros(Nx, Ny, Nz, param.nimfs);
Residue.u = u; Residue.v = v;

H = zeros(Nx,Ny,Nz,2);
H1 = zeros(Nx,Ny,Nz,2);
mse = zeros(2,1);

Windows = zeros(7,param.nimfs);

sift_cnt = zeros(1,param.nimfs);
imf = 1;
stopflag = 1;
    while(imf <= param.nimfs && stopflag)
        %Initialising intermediary IMFs
        H(:,:,:,1) = Residue.u; H(:,:,:,2) = Residue.v;

        sift_stop = 0; %flag to control sifting loop
        
        Combined = H(:,:,:,1)/sqrt(2) + H(:,:,:,2)/sqrt(2); %Combining two signals with equal weights        
        [Maxima,MaxPos,Minima,MinPos] = MinimaMaxima3D(Combined,1,1,[],[]);  %Obtaining extrema of combined signal
        
        %Checking whether there are too few extrema in the IMF
        if (nnz(Maxima) < 3 || nnz(Minima) < 3)
            warning('Fewer than three extrema found in extrema map. Stopping now...');
            break;
        end
        
        %Window size determination by delaunay triangulation
        Windows(:,imf) = filter_size3D(MaxPos,MinPos,param.type);        
        w_sz = Windows(param.type,imf); %extracting window size chosen by input parameter
        
        if~(w_sz)
           warning('EMD3D3V has stopped because the Delaunay Triangulation could not be created (collinear points)'); 
           param.nimfs = imf-1;
           break;
        end
        
        %Begin sifting iteration
        while~(sift_stop)            
            sift_cnt(imf) = sift_cnt(imf) + 1; %Incrementing sift counter
            
            %Entering parallel sift calculations
            
            parfor i=1:2
               H1(:,:,:,i) = Sift(H(:,:,:,i),w_sz);
               
               mse(i) = immse(H1(:,:,:,i),H(:,:,:,i));
            end
                       
            %Stop condition checks       
            if (mse(1)<param.tol && mse(2)<param.tol && sift_cnt(imf)~=1)
                sift_stop = 1;
            end
            
            H = H1;         
        end
        
        %Storing IMFs
        IMF.u(:,:,:,imf) = H(:,:,:,1); IMF.v(:,:,:,imf) = H(:,:,:,2);

        %Subtracting from Residual Signals
        Residue.u = Residue.u - IMF.u(:,:,:,imf);
        Residue.v = Residue.v - IMF.v(:,:,:,imf);
        
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
    [Results.IO.u,Results.Error.u] = Orth_index3D(u,IMF.u,Residue.u);
    [Results.IO.v,Results.Error.v] = Orth_index3D(v,IMF.v,Residue.v);
    
    switch(param.plot)
        case 'on'
            Plot_results(u,v,Results,param)
    end
end

function H1 = Sift(H,w_sz)

%Envelope Generation
[Env_max,Env_min] = OSF(H,w_sz);

%padding
Env_med = Pad_smooth(Env_max,Env_min,w_sz);

%Subtracting from residue
H1 = H - Env_med;
                
end

function [Max,Min] = OSF(H,w_sz)
%Order statistics filtering to determine maximum and minmum envelopes
            Max = Separable_ordfilt3(H, 'max', w_sz); %Max envelope
            Min = Separable_ordfilt3(H, 'min', w_sz); %Min envelope 
end

function Env_med = Pad_smooth(Env_max,Env_min,w_sz)
h = floor(w_sz/2);

%Padding
temp = padarray(Env_max,[h h],'replicate');
temp1 = permute(temp,[3 2 1]); %interchanging dimensions
temp = padarray(temp1,[h 0],'replicate');
Env_maxp = permute(temp,[3 2 1]); %restoring dimensions

temp = padarray(Env_min,[h h],'replicate');
temp1 = permute(temp,[3 2 1]); %interchanging dimensions
temp = padarray(temp1,[h 0],'replicate');
Env_minp = permute(temp,[3 2 1]); %restoring dimensions

%Smoothing

temp1 = movmean(Env_maxp,w_sz,3,'endpoints','discard');
temp2 = movmean(temp1,w_sz,2,'endpoints','discard');
Env_maxs = movmean(temp2,w_sz,1,'endpoints','discard');

temp1 = movmean(Env_minp,w_sz,3,'endpoints','discard');
temp2 = movmean(temp1,w_sz,2,'endpoints','discard');
Env_mins = movmean(temp2,w_sz,1,'endpoints','discard');

%Calculating mean envelope
Env_med = (Env_maxs + Env_mins)./2;

end

function Plot_results(u,v,Results,param)
% default plot attributes
set(groot,'defaultaxesfontname','times');
set(groot,'defaultaxesfontsize',12);
set(groot,'defaulttextInterpreter','latex');
set(groot,'defaultLineLineWidth',2);

Colour = parula;
nslice  = param.nslice;

figure(1)   
        subplot(1,2,1)
        TIMF_plot(u,Colour,nslice,0,'Signal','u');
        subplot(1,2,2)
        TIMF_plot(v,Colour,nslice,0,'Signal','v');



    for i=1:param.nimfs
     figure(i+1)   
        subplot(1,2,1)
        TIMF_plot(Results.IMF.u(:,:,:,i),Colour,nslice,i,'IMF','u');
        subplot(1,2,2)
        TIMF_plot(Results.IMF.v(:,:,:,i),Colour,nslice,i,'IMF','v');
    end
    
    figure(i+2)
    subplot(1,2,1)
        TIMF_plot(Results.Residue.u,Colour,nslice,0,'Residue','u');
        subplot(1,2,2)
        TIMF_plot(Results.Residue.v,Colour,nslice,0,'Residue','v');
end