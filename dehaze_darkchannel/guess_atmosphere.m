function [ A ] = guess_atmosphere( I, D )

[h, w, ~] = size(I);
N = h * w;

Nsearch = floor(N * 0.01);

I = reshape(I, N, 3);

[~, ind] = sort(D(:), 'descend');

A = mean(I(ind(1:Nsearch), :), 1);

end

