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
sliceViewer(im2double(abs(Vy)));
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
s = 1;
m = (n+1)/2;

G = @(x,m,s) 1/(s*sqrt(2*pi)) .* exp(-(x-m).^2 ./ (2*s^2));
dGdx = @(x,m,s) -(x-m)/s^2 .* G(x,m,s);

g_filter = zeros(1,n,1);
g_filter(1,:,1) = G(1:n,m,s);
dg_filter(1,:,1) = dGdx(1:n,m,s);

% x gauss
I_filterX = imfilter(I,permute(g_filter,[1,3,2]));
I_filterX = imfilter(I_filterX,permute(g_filter,[2,1,3]));
%sliceViewer(I_filter);
I_filterX = imfilter(I_filterX,dg_filter);

% y gauss
g_filter = zeros(1,n,1);
g_filter(1,:,1) = G(1:n,m,s);
dg_filter(1,:,1) = dGdx(1:n,m,s);

I_filterY = imfilter(I,g_filter);
I_filterY = imfilter(I_filterY,permute(g_filter,[1,3,2]));
I_filterY = imfilter(I_filterY,permute(dg_filter,[2,1,3]));


% t gauss
I_filterT = imfilter(I,g_filter);
I_filterT = imfilter(I_filterT,permute(g_filter,[2,1,3]));
I_filterT = imfilter(I_filterT,permute(dg_filter,[1,3,2]));

%%
N = 9;
A = zeros(N^2,2);
b = zeros(N^2,1);
%[x,y,t] = deal(130,100,32);
t = 56;
% x = 5:10:250;
% y = 5:10:250;
x = 10:5:240;
y = 10:5:240;
count = 1;
for i = x
    for j = y
        
        idx_x = (1:N)-(N+1)/2 + i;
        idx_y = (1:N)-(N+1)/2 + j;
        A(:,1) = reshape(I_filterX(idx_x,idx_y,t),[N^2,1]);
        A(:,2) = reshape(I_filterY(idx_x,idx_y,t),[N^2,1]);
        b(:) = reshape(I_filterT(idx_x,idx_y,t),[N^2,1]);
        tmp = A\(-b);
        if norm(tmp) < 1
            continue
        end
        dX(count) = tmp(1);
        dY(count) = tmp(2);
        X(count) = i;
        Y(count) = j;
        count = count + 1;
    end
end


figure(11)
imshow(I(:,:,t))
hold on
alpha = 1;
quiver(Y,X,alpha.*dX,alpha.*dY,2,LineWidth=2)
%%
figure(5)
sliceViewer(abs(I_filterX));
figure(6)
sliceViewer(abs(I_filterY));
figure(8)
sliceViewer(abs(I_filterT));

%%
% hzm = zeros(1,1,3);
% hzm(1,1,:) = [1,0,-1];
% hy = zeros(3,1,1);
% hy(:,1,1) = [1,2,1];
% hx = zeros(1,3,1);
% hx(1,:,1) = [1,2,1];
% 
% hx'*hy
