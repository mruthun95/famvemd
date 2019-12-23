% Purpose: 
% -To perform EMD on n(3-16) channels of 1 dimensional data
%
% Input: 
% - u: Signal of size (length, number of channels)
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

function Results = EMD1DNV(u,t,param) 

%Reading signal characteristics
[Nx] = size(u,1); %Signal dimensions
    
 n = size(u,2); %number of channels

if(~all(ismember(param.type,[1,2,3,4,5,6,7])))
    error('Please enter a valid window size type')
end
 
if ~isfield(param,'nimfs')
    param.nimfs = 10;
end

if ~isfield(param,'tol')
    C = min(rms(u,1));
    param.tol = C*0.001; % 0.1% of the minimum signal amplitude
end

if ~isfield(param,'type')
    param.type = 5;
end

if ~isfield(param,'plot')
    param.plot = 'off';
end

%Initialisations
IMF     = zeros(Nx,n,param.nimfs); 
H1      = zeros(Nx,n);
mse     = zeros(n,1);

Windows = zeros(7,param.nimfs);

sift_cnt = zeros(1,param.nimfs);
imf = 1;

Residue = u;

    while(imf <= param.nimfs)
        %Initialising intermediary IMFs
        H = Residue;

        sift_stop = 0; %flag to control sifting loop
            
        Combined = sum(H/sqrt(n),2); %Combining two signals with equal weights 
        
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
            
            %Entering parallel sift calculations
            
            for i=1:n
               H1(:,i) = Sift(H(:,i),w_sz);
               
               mse(i) = immse(H1(:,i),H(:,i));
            end
            
            %Stop condition checks       
            if ( all(mse<param.tol) && sift_cnt(imf)~=1)
                sift_stop = 1;
            end
            
            H = H1 ;    
        end
        
        %Storing IMFs
        IMF(:,:,imf) = H;

        %Subtracting from Residual Signals
        Residue = Residue - IMF(:,:,imf);
       
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
    [Results.IO,Results.Error] = Orth_index(u,IMF,Residue);

    switch(param.plot)
        case 'on'
            Plot_results(u,t,Results)
    end
end

function Windows = filter_size1D(imax, imin, type)
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

edge_len_max = diff(sort(imax));
edge_len_min = diff(sort(imin));
    
    
        %Window size calculations
        
        d1 = min( min(edge_len_max) , min(edge_len_min) );
        d2 = max( min(edge_len_max) , min(edge_len_min) );
        d3 = min( max(edge_len_max) , max(edge_len_min) );
        d4 = max( max(edge_len_max) , max(edge_len_min) );
        d5 = (d1+d2+d3+d4)/4 ;
        d6 = median([edge_len_min; edge_len_max]);
        d7 = mode([edge_len_min; edge_len_max]);
        
        Windows = [d1, d2, d3, d4, d5, d6, d7];

%making sure w_size is an odd integer
Windows = 2*(floor(Windows./2))+1;
         
if(Windows(type)<3)
    warning('WARNING: Calculated Window size less than 3');
    warning('Overriding calculated value and setting window size = 3');
    Windows(type) = 3;
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
            Max = Ordfilt1(H, 'max', w_sz); %Max envelope
            Min = Ordfilt1(H, 'min', w_sz); %Min envelope 
                
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

function Env_med = Pad_smooth(Env_max,Env_min,w_sz)
h = floor(w_sz/2);

%Padding
%u
Env_maxp = padarray(Env_max,[h 0],'symmetric');
Env_minp = padarray(Env_min,[h 0],'symmetric');


%Smoothing
%u
Env_maxs = movmean(Env_maxp,w_sz,1,'endpoints','discard');
Env_mins = movmean(Env_minp,w_sz,1,'endpoints','discard');

%Calculating mean envelope
Env_med = (Env_maxs + Env_mins)./2;
end

function [IO,Error] = Orth_index(Signal,IMF,Residue)
% Purpose: 
% To calculate the index of orthogonality of a decomposition and its mean
% squared error

n = size(Signal,2); %Number of channels
n_imf = size(IMF,3);
numerator = zeros(size(Signal));
I = sum(IMF,3) + Residue;

Error.global = zeros(1,n);

for i=1:n
Error.global(i) = immse(I(:,i),Signal(:,i));
end


for j = 1:n_imf
    for k = 1:n_imf
        if(j~=k)
           numerator = numerator + IMF(:,:,j).*IMF(:,:,k);
        end
    end
end

IO.map = numerator/(sum(Signal.^2)); %wrong
IO.global = sum(IO.map);
end

function Plot_results(u,t,Results)
% default plot attributes
set(groot,'defaultaxesfontname','times');
set(groot,'defaultaxesfontsize',12);
set(groot,'defaulttextInterpreter','tex');
set(groot,'defaultLineLineWidth',2);

n = size(u,2); %Number of channels
n_imfs = size(Results.IMF,3);

figure(1)   
skip = 0;
    for j = 1:n
       for i=1:n_imfs+2 
          
        if i==1 
            strng = sprintf('%d',j);
            subplot(n,n_imfs+2,i+skip)
            IMF_plot(u(:,j),t,0,'Channel',strng);
        elseif i>1 && i<n_imfs+2
            strng = sprintf('Channel %d',j);
            subplot(n,n_imfs+2,i+skip)
            IMF_plot(Results.IMF(:,j,i-1),t,i-1,'IMF',strng);
        else
            strng = sprintf('Channel %d',j);
            subplot(n,n_imfs+2,i+skip)
            IMF_plot(Results.Residue(:,j),t,0,'Residue',strng);
        end
       end
       skip = skip + n_imfs+2;
    end
    
end

function IMF_plot(signal,t,imf,name1,name2)    

    plot(t,signal,'-k');
    xlabel('t');
    set(gca,'TickLabelInterpreter','tex')
    switch(name1)
        case 'IMF'
            title(sprintf('%s %d %s',name1,imf,name2));
        case 'Channel'
            title(sprintf('%s %s',name1,name2));
        case 'Residue'
            title(sprintf('%s %s',name1,name2));
    end
    
    
end

