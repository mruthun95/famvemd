function [IO,Error] = Orth_index2D(Signal,IMF,Residue)
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