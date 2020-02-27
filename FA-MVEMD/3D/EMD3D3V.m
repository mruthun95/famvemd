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
        Windows(:,imf) = filter_size3D(MaxPos,MinPos,param.type);        
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
    [Results.IO.u,Results.Error.u] = Orth_index3D(u,IMF.u,Residue.u);
    [Results.IO.v,Results.Error.v] = Orth_index3D(v,IMF.v,Residue.v);
    [Results.IO.w,Results.Error.w] = Orth_index3D(w,IMF.w,Residue.w);
    
    switch(param.plot)
        case 'on'
            Plot_results(u,v,w,Results,param)
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