function [ Q ] = guided_filter( I, P, psz )

[h, w, ch] = size(I);
N = h * w;
pw = floor(psz/2);

Pmean = colfilt(padarray(P, [pw pw], nan), [psz psz], 'sliding', @nanmean);
Pmean = Pmean(1+pw:end-pw, 1+pw:end-pw);

A = zeros(size(I));
B = zeros(h, w);

% naive implementation; need tests
for k = 1:N
	% k-th window
	[r, c] = ind2sub([h, w], k);
	rk = max(1, r-pw):min(h, r+pw);
	ck = max(1, c-pw):min(w, c+pw);
	
	Nk = numel(rk) * numel(ck);
	Ik = reshape(I(rk, ck, :), [Nk ch])';
	Pk = reshape(P(rk, ck), [1 Nk]);
	mu = mean(Ik, 2);
	S = Ik * Ik' / Nk - mu * mu';
	ak = (S + eps * eye(ch)) \ (sum(Ik .* repmat(Pk, 3, 1), 2) / Nk - mu * Pmean(k));
	bk = Pmean(k) - ak' * mu;
	
	A(r, c, :) = ak;
	B(r, c) = bk;
end

Amean = zeros(size(I));
Bmean = colfilt(padarray(B, [pw pw], nan), [psz psz], 'sliding', @nanmean);
Bmean = Bmean(1+pw:end-pw, 1+pw:end-pw);

for k = 1:N
	[r, c] = ind2sub([h, w], k);
	rk = max(1, r-pw):min(h, r+pw);
	ck = max(1, c-pw):min(w, c+pw);
	
	Nk = numel(rk) * numel(ck);
	
	Amean(r, c, :) = mean(reshape(I(rk, ck, :), [Nk 3]));
end

Q = sum(Amean .* I, 3) + Bmean;

end

