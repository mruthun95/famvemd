function [IO,Error] = Orth_index3D(Signal,IMF,Residue)
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