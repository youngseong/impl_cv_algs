function [ L ] = construct_matting_laplacian( I, wsz )

if ~exist('wsz', 'var')
	wsz = 3;
end

[h, w, ~] = size(I);
N = h * w;
pw = floor(wsz / 2);
Np = wsz^4;

ii = zeros(Np * N, 1);
jj = zeros(Np * N, 1);
vv = zeros(Np * N, 1);
len = 0;

epsilon = 1e-4;

for k = 1:N
	[r, c] = ind2sub([h, w], k);
	% k-th window
	rk = max(1, r-pw):min(h, r+pw);
	ck = max(1, c-pw):min(w, c+pw);
	Nk = numel(rk) * numel(ck);
	Iwk = reshape(I(rk, ck, :), [Nk 3])';
	mu = mean(Iwk, 2);
	C = Iwk * Iwk' / Nk - mu * mu';
	Idiff = Iwk - repmat(mu, [1 Nk]);
	
	vals = eye(Nk) - 1 / Nk * (1 + Idiff' / (C + epsilon / Nk * eye(3)) * Idiff);
	
	[Xk, Yk] = meshgrid(ck, rk);
	inds = sub2ind([h w], Yk(:), Xk(:));
	
	subrange = len + (1:Nk^2);
	
	ii(subrange) = repmat(inds, 1, Nk);
	jj(subrange) = repmat(inds', Nk, 1);
	vv(subrange) = vals(:);
	len = len + Nk^2;
end

L = sparse(ii(1:len), jj(1:len), vv(1:len));

end

