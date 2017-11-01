function [ J ] = recover_radiance( I, A, t, t0 )

[h, w, ~] = size(I);

A = repmat(reshape(A, [1 1 3]), h, w);
J = (I - A) ./ repmat(max(t, t0), 1, 1, 3) + A;

end

