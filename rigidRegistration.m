% implementation of two 2D point-based non-rigid registration methods
% 1) second order polynomials
% 2) thin-plate splines


close all
clear, clc
 
load('HW7_practice5.mat');
%file = 'Please input the mat:\n';
%path = input(file,'s');
%load(path);
 
%show original Image
win = figure;
figure(win);
subplot(2,2,1);
imshow(source, []);
hold on;
plot(xc1, yc1, 'r*');
hold off;
title('Source Image');
 
figure(win);
subplot(2,2,2);
imshow(target, []); 
hold on;
plot(xc2, yc2, 'g*');
hold off;
title('Target Image');
 
%method 1:
if method == 1
    
    sourceImage = [xc1(:), yc1(:)];
    targetImage = [ones(size(xc2(:))), xc2(:), yc2(:), xc2(:).*yc2(:), xc2(:).^2, yc2(:).^2];
    transformMatrix = targetImage \ sourceImage;
    %disp(transformMatrix);
    %transformMatrix = pinv(targetImage) * sourceImage; 
    disp(transformMatrix)
    
    [X, Y] = meshgrid(1:5:256, 1:5:256);
    grid=zeros(size(sourceImage));
    grid(X,1:256) = 1;
    grid(1:256, Y) = 1;
 
    [xx, yy] = meshgrid(1:256, 1:256);
    xx = xx(:);
    yy = yy(:);
 
    template = [ones(size(xx)), xx, yy, xx.*yy, xx.^2, yy.^2] * transformMatrix; 
    xSource = template(:, 1);
    ySource = template(:, 2); 
 
    x = reshape(xSource, size(source));
    y = reshape(ySource, size(source));
 
    regImage = interp2(source, x,y,'bilinear',0);
    gridField = interp2(grid, x, y, 'bilinear',0);
    
    figure(win);
    subplot(2,2,3);
    imshow(regImage, []);
    hold on;
    plot(xc2, yc2, 'g*');
    hold off;
    title('registrated Image');
 
    figure(win);
    subplot(2,2,4);
    imshow(gridField, []); 
    title('Deformable Field');    
end
 
if method == 2
    
    sourceImage = zeros(size(xc1,2)+3, 2);
    sourceImage(1:size(xc1,2), 1) = xc1(:);
    sourceImage(1:size(xc1,2), 2) = yc1(:);
    targetImage = [xc2(:), yc2(:)];
    dis = pdist2(targetImage,targetImage,'euclidean').^2;
    dis(dis==0) = 1;
    r = dis .* log(dis);
    
    A = zeros(size(dis,1),size(xc1,2)+3);
    A(1:size(xc1,2),1) = 1;
    A(1:size(xc1,2),2) = xc2(:);
    A(1:size(xc1,2),3) = yc2(:);
    A(:,4:end) = r;
    A(size(xc1,2)+1,4:end) = 1;
    A(size(xc1,2)+2,4:end) = xc2(:);
    A(size(xc1,2)+3,4:end) = yc2(:); 
    
    %transformMatrix = A \ sourceImage;
    %disp(transformMatrix);
    transformMatrix = pinv(A) * sourceImage; 
    disp(transformMatrix)
    
    [X, Y] = meshgrid(1:5:256, 1:5:256);
    grid=zeros(size(sourceImage));
    grid(X,1:256) = 1;
    grid(1:256, Y) = 1;
 
    [xx, yy] = meshgrid(1:256, 1:256);
    xx = xx(:);
    yy = yy(:);
    meshPoint = [xx,yy];
    
    %dis = pdist2(meshPoint,targetImage,'euclidean').^2;
    sourece = [xc2(:), yc2(:)];
    dis = pdist2(meshPoint,sourece,'euclidean').^2;
    r = dis .* log(dis);
 
    meshA = zeros(size(dis,1),size(xc1,2)+3);
    meshA(:,1) = 1;
    meshA(:,2) = xx;
    meshA(:,3) = yy;
    meshA(:,4:end) = r;
 
    template = meshA * transformMatrix;
    xSource = template(:, 1);
    ySource = template(:, 2);
    x = reshape(xSource, size(source));
    y = reshape(ySource, size(source));
    regImage = interp2(source, x,y,'bilinear',0);
    gridField = interp2(grid, x, y, 'bilinear',0);
      
    figure(win);
    subplot(2,2,3);
    imshow(regImage, []);
    hold on;
    plot(xc2, yc2, 'g*');
    hold off;
    title('registrated Image');
 
    figure(win);
    subplot(2,2,4);
    imshow(gridField, []); 
    title('Deformable Field');
    
end

