%%
clear, close all
%%
clear
%%
I = zeros(256,144,82-40+1);
A = imresize( im2double(rgb2gray(imread('video/00005.png'))), 1/5 );
for k = 40:82
    I(:,:,k) = imresize( im2double(rgb2gray(imread(sprintf('video/%05d.png',k)))), 1/5 );
end

figure(1)
sliceViewer(I);

%%
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
I_filterT = imfilter(I,g_filter,'replicate');
I_filterT = imfilter(I_filterT,permute(g_filter,[2,1,3]),'replicate');
I_filterT = imfilter(I_filterT,permute(dg_filter,[1,3,2]),'replicate');

%%
for t = 45:77
    N = 9;
    A = zeros(N^2,2);
    b = zeros(N^2,1);
    %[x,y,t] = deal(130,100,32);
    % t = 56;
    % x = 5:10:250;
    % y = 5:10:250;
    x = 10:5:240;
    y = 10:5:135;
    count = 1;
    for i = x
        for j = y
            
            idx_x = (1:N)-(N+1)/2 + i;
            idx_y = (1:N)-(N+1)/2 + j;
            A(:,1) = reshape(I_filterX(idx_x,idx_y,t),[N^2,1]);
            A(:,2) = reshape(I_filterY(idx_x,idx_y,t),[N^2,1]);
            b(:) = reshape(I_filterT(idx_x,idx_y,t),[N^2,1]);
            tmp = A\(-b);
            if norm(tmp) < 0.2
                continue
            end
            dX(count) = tmp(1);
            dY(count) = tmp(2);
            X(count) = i;
            Y(count) = j;
            count = count + 1;
        end
    end
    
    
    figure()
    imshow(I(:,:,t))
    hold on
    alpha = 1;
    quiver(Y,X,alpha.*dX,alpha.*dY,2,LineWidth=2)
    H = gcf;
    H.WindowState = 'maximized';
%     exportgraphics(gcf,sprintf('visualization/frame_%02d.png',t),'Resolution',1000)
    saveas(gcf, sprintf('video_analysis/frame_%02d.png',t));
end

%%

for k = 45:77
    filename = sprintf('video_analysis/frame_%02d.png',k);
    vis(:,:,k) = im2double(rgb2gray(imread(filename)));
end
sliceViewer(vis);
