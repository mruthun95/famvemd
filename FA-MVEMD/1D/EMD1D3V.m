% Purpose: 
% -To perform EMD on 3 channels of 1 dimensional data
%
% Input: 
% - u: Signal 1
% - v: Signal 2
% - w: Signal 3
% - param
%   -nimfs: Number of IMFs to be extracted 
%   -tol: Sifting tolerance value
%   -type: type of window size to be used
%   -plot: 'on' to plot results, default hides IMF plots
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

function Results = EMD1D3V(u,v,w,t,param) 

%Reading signal characteristics
[Nx] = length(u); %Signal dimensions
B    = length(v); %Signal dimensions
C    = length(w); %Signal dimensions

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

if(~all(Nx==B) ||  ~all(Nx==C))
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
IMF.u = zeros(Nx, param.nimfs); 
IMF.v = zeros(Nx, param.nimfs);
IMF.w = zeros(Nx, param.nimfs);
Residue.u = u; Residue.v = v; Residue.w = w;

Windows = zeros(7,param.nimfs);

sift_cnt = zeros(1,param.nimfs);
imf = 1;

    while(imf <= param.nimfs)
        %Initialising intermediary IMFs
        H.u = Residue.u; H.v = Residue.v; H.w = Residue.w;

        sift_stop = 0; %flag to control sifting loop
        
        Combined = H.u/sqrt(3) + H.v/sqrt(3) + H.w/sqrt(3) ; %Combining two signals with equal weights        
        [Maxima,MaxPos,Minima,MinPos] = extrema(Combined);  %Obtaining extrema of combined signal
        
        
        %Checking whether there are too few extrema in the IMF
        if (nnz(Maxima) < 3 || nnz(Minima) < 3)
            warning('Fewer than three extrema found in extrema map. Stopping now...');
            break;
        end
        
        %Window size determination by delaunay triangulation
        Windows(:,imf) = filter_size1D(MaxPos,MinPos,param.type);        
        w_sz = Windows(param.type,imf); %extracting window size chosen by input parameter
        
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
            
%             osfplot(t, H, Env);
            
            %Subtracting from residue
            H1.u = H.u - Env.u.med; H1.v = H.v - Env.v.med; H1.w = H.w - Env.w.med;     
            
%             projenv_IMF_plot(t, Combined, MaxPos,MinPos)
             
            %Stop condition checks
            mse_u = immse(H1.u,H.u); mse_v = immse(H1.v,H.v); mse_w = immse(H1.w,H.w);       
            if (mse_u<param.tol && mse_v<param.tol && mse_w<param.tol && sift_cnt(imf)~=1)
                sift_stop = 1;
            end
            
            H.u = H1.u; H.v = H1.v; H.w = H1.w;           
        end
        
        %Storing IMFs
        IMF.u(:,imf) = H.u; IMF.v(:,imf) = H.v; IMF.w(:,imf) = H.w;

        %Subtracting from Residual Signals
        Residue.u = Residue.u - IMF.u(:,imf);
        Residue.v = Residue.v - IMF.v(:,imf);
        Residue.w = Residue.w - IMF.w(:,imf);
        
        %Incrementing IMF counter
        imf = imf + 1;
        
    end
    
     %Checking for oversifting
    if(any(sift_cnt>=5*ones(1,param.nimfs)))
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
    [Results.IO.w,Results.Error.w] = Orth_index(w,IMF.w,Residue.w);

    switch(param.plot)
        case 'on'
            Plot_results(u,v,w,t,Results,param)
    end
end

function Env = OSF(H,w_sz)
%Order statistics filtering to determine maximum and minmum envelopes
            Env.u.max = Ordfilt1(H.u, 'max', w_sz); %Max envelope u
            Env.u.min = Ordfilt1(H.u, 'min', w_sz); %Min envelope u
             
            Env.v.max = Ordfilt1(H.v, 'max', w_sz); %Max envelope v
            Env.v.min = Ordfilt1(H.v, 'min', w_sz); %Min envelope v
            
            Env.w.max = Ordfilt1(H.w, 'max', w_sz); %Max envelope w
            Env.w.min = Ordfilt1(H.w, 'min', w_sz); %Min envelope w
end

function Env = Pad_smooth(Env,w_sz)
h = floor(w_sz/2);

%Padding
%u
Env.u.maxp = padarray(Env.u.max,[h 0],'symmetric');
Env.u.minp = padarray(Env.u.min,[h 0],'symmetric');
%v
Env.v.maxp = padarray(Env.v.max,[h 0],'symmetric');
Env.v.minp = padarray(Env.v.min,[h 0],'symmetric');
%w
Env.w.maxp = padarray(Env.w.max,[h 0],'symmetric');
Env.w.minp = padarray(Env.w.min,[h 0],'symmetric');

%Smoothing
%u
Env.u.maxs = movmean(Env.u.maxp,w_sz,1,'endpoints','discard');
Env.u.mins = movmean(Env.u.minp,w_sz,1,'endpoints','discard');
%v
Env.v.maxs = movmean(Env.v.maxp,w_sz,1,'endpoints','discard');
Env.v.mins = movmean(Env.v.minp,w_sz,1,'endpoints','discard');
%w
Env.w.maxs = movmean(Env.w.maxp,w_sz,1,'endpoints','discard');
Env.w.mins = movmean(Env.w.minp,w_sz,1,'endpoints','discard');

end

function [IO,Error] = Orth_index(Signal,IMF,Residue)
% Purpose: 
% To calculate the index of orthogonality of a decomposition and its mean
% squared error

n_imf = size(IMF,2);
numerator = zeros(size(Signal));
I = sum(IMF,2) + Residue;

Error.map = (Signal-I)./Signal;
Error.global = immse(I,Signal);

for j = 1:n_imf
    for k = 1:n_imf
        if(j~=k)
           numerator = numerator + IMF(:,j).*IMF(:,k);
        end
    end
end

IO.map = numerator/(sum(Signal.^2)); %wrong
IO.global = sum(IO.map);
end

function Plot_results(u,v,w,t,Results,param)
% default plot attributes
set(groot,'defaultaxesfontname','times');
set(groot,'defaultaxesfontsize',12);
set(groot,'defaulttextInterpreter','tex');
set(groot,'defaultLineLineWidth',2);

figure(1)   
        subplot(3,1,1)
        IMF_plot(u,t,0,'Channel','1');
        subplot(3,1,2)
        IMF_plot(v,t,0,'Channel','2');
        subplot(3,1,3)
        IMF_plot(w,t,0,'Channel','2');


     for i=1:param.nimfs   
     figure(i+1)    
        subplot(3,1,1)
        IMF_plot(Results.IMF.u(:,i),t,i,'IMF','Channel 1');
        subplot(3,1,2)
        IMF_plot(Results.IMF.v(:,i),t,i,'IMF','Channel 2');
        subplot(3,1,3)
        IMF_plot(Results.IMF.w(:,i),t,i,'IMF','Channel 3');
     end
    
    figure(i+2)  
    subplot(3,1,1)
    IMF_plot(Results.Residue.u,t,0,'Residue','Channel 1');
    subplot(3,1,2)
    IMF_plot(Results.Residue.v,t,0,'Residue','Channel 2');
    subplot(3,1,3)
    IMF_plot(Results.Residue.w,t,0,'Residue','Channel 3');
end