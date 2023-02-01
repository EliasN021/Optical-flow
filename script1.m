%%
clear, close all
%%
I = zeros(256,256,64);
for k = 1:64
    filename = sprintf('toyProblem_F22/frame_%02d.png',k);
    I(:,:,k) = im2double(rgb2gray(imread(filename)));
end

figure(1)
sliceViewer(I);


Vx = I(2:end,:,:) - I(1:end-1,:,:);
Vy = I(:,2:end,:) - I(:,1:end-1,:);
Vt = I(:,:,2:end) - I(:,:,1:end-1);

figure(2)
sliceViewer(im2double(abs(Vx)));
figure(7)
%%
PWx = [1,0,-1; 1,0,-1; 1,0,-1];
PWy = [1,1,1; 0,0,0; -1,-1,-1];
I_PWx = imfilter(I,PWx);
I_PWy = imfilter(I,PWy);

figure(3)
sliceViewer(im2double(abs(I_PWx)));
figure(4)
sliceViewer(im2double(abs(I_PWy)));

%%
% d/dx(G(x)*G(y)) = G'(x)*G(y) + G(x)*G'(y)
n = 5;
s = 4;
m = (n+1)/2;

G = @(x,m,s) 1/(s*sqrt(2*pi)) .* exp(-(x-m).^2 ./ (2*s^2));
dGdx = @(x,m,s) -(x-m)/s^2 .* G(x,m,s);

g_filter = zeros(1,n,1);
g_filter(1,:,1) = G(1:n,m,s);
dg_filter(1,:,1) = dGdx(1:n,m,s);

%% x gauss
I_filter = imfilter(I,permute(g_filter,[1,3,2]));
I_filter = imfilter(I_filter,permute(g_filter,[2,1,3]));
%sliceViewer(I_filter);
I_filter = imfilter(I_filter,dg_filter);

figure(5)
sliceViewer(abs(I_filter));
