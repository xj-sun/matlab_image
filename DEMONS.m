
% implement a 2D version of Thirion’s DEMONS algorithm

close all
clear, clc
load('TestHw8_I6L3_15.mat');

figure;
win = [];
win(1) = subplot(2,2,1); 
win(2) = subplot(2,2,2); 
win(3) = subplot(2,2,3); 
win(4) = subplot(2,2,4); 
colormap('gray')

moveImage = double(dynamic);
fixImage = double(static);
transformX = zeros(size(dynamic));
transformY = zeros(size(dynamic));

alpha = 2;
kernelSize = 60;
iter = 300;
for level = numlevel : -1 : 1
    
    size = 256 / 2^(level-1);
    mm = imresize(moveImage, [size, size]);
    ff = imresize(fixImage, [size, size]);
    transformX = imresize(transformX, [size, size]);
    transformY = imresize(transformY, [size, size]);
    
    [X, Y] = meshgrid(1:size);
    [Gx, Gy] = gradient(ff);
    
    imagesc(mm, 'Parent', win(1)); title(win(1), 'downsample dynamic'); axis(win(1), 'square');
    imagesc(ff, 'Parent', win(2)); title(win(2), 'downsample static'); axis(win(2), 'square');
    
    Skernel = fspecial('gaussian', [kernelSize kernelSize], 10);
    for itt = 1 : iter
        
        x = X - transformX;
        y = Y - transformY;
        mmCur = interp2(X, Y, mm, x, y);
        disAb = abs(mmCur - ff);
        imagesc(mmCur,'Parent',win(3)); title(win(3), 'registered process'); axis(win(3), 'square')
        imagesc(disAb,'Parent',win(4)); title(win(4), 'absolutr distance'); axis(win(4), 'square');
        pause(0.01);   
             
        distance = ff - mmCur;  % Demon force
        Ux = -(distance.*Gx)./ ((Gx.^2+Gy.^2) + distance.^2);
        Uy = -(distance.*Gy)./ ((Gx.^2+Gy.^2) + distance.^2);
        
        Ux(isnan(Ux)) = 0; 
        Uy(isnan(Uy)) = 0; 
        Ux = imfilter(Ux, Skernel);
        Uy = imfilter(Uy, Skernel);
        transformX = transformX + Ux;
        transformY = transformY + Uy;
    end
    
    transformX = transformX .* 2; 
    transformY = transformY .* 2;      
end
