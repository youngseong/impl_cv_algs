% input: color image (N*M*3)
% output: J(N*M*3), A(3), t(N*M)

close all; clear; clc;

I = double(imread('tiananmen1.png'))/255;
% I = imresize(I, .25);

[h, w, ~] = size(I);
N = h * w;
psz = 15;
om = .95;
eps = 1e-7;
lambda = 1e-4;

D = compute_dark_channel(I, psz);
A = guess_atmosphere(I, D);
% A = max(reshape(I, [h*w, 3]));

% (12)
t_ = 1 - om * compute_dark_channel(I ./ repmat(reshape(A, [1 1 3]), h, w), psz);

tic;
L = construct_matting_laplacian(I);
toc;

% (17)
t = (L + lambda * speye(N)) \ (lambda * t_(:));
t = reshape(t, [h w]);
figure, imagesc(t);

% find J
J = recover_radiance(I, A, t, .1);

figure, imshow(I);
figure, imshow(J);

imwrite(I, '_in.png')
imwrite(J, '_out.png')