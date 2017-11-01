function [ D ] = compute_dark_channel( I, psz )

pw = floor(psz/2);

D = min(I, [], 3);
D = padarray(D, [pw pw], nan);
D = colfilt(D, [psz psz], 'sliding', @nanmin);
D = D(pw+1:end-pw, pw+1:end-pw);

end

