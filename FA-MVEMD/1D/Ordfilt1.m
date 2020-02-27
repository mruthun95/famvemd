
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