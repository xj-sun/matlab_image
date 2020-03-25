
% implement the automated model-based bias field correction for MR images described
%in the paper by Van Leemput et al

close all; clear;

file=input('What is the name of your file?', 's')
%file='testhw3_1';
IStruct=load(file);
img = IStruct.ima_att;
figure(1);
subplot(2,2,1);
imagesc(img); 
title('Original image');

pixVal = img(:);

[counts, centers]=hist(pixVal, 0:255);
f1=counts/sum(counts);
subplot(2,2,2);
hold on
plot(f1);
title('Original image');

x = 1:size(img,1);
y = 1:size(img,2);
[X,Y] = meshgrid(x,y);
subplot(2,2,3);
mesh(X,Y,img);
grid on
title('Original image');


pause(3);

maxiter=40;


[M,N,P] = size(img);
biasParams = [0,0,0,0,0,0];
biasSurface = biasSurf(img, biasParams);
figure(2)

subplot(3,2,1);plot(1:M, img(:,128),'b-');title('Original image (:,128)');
subplot(3,2,3);plot(1:M, img(128,:),'b-');title('Original image (128,:)');

pause(4)
reImage = reshape(img, [P, M*N]);


model= init_values(IStruct);
label = E_step(reImage, model);
correct = img - biasSurface;
[X, Y] = meshgrid(1:M, 1:N);

for iter = 1: maxiter
     
    model = M_step(correct,label);
    
    
    biasParams = getC(model, label, reImage,img);
    biasSurface = biasSurf(img, biasParams);
    correct = img - biasSurface;
    
    reCorrect=correct(:);
    [counts, centers]=hist(reCorrect, 0:255);
    f2=counts/sum(counts);
    


    figure(3);
    subplot(2,2,1);plot(f1,'b');hold on; plot(f2,'r');title('Original and estimated distibution');
    hold off;
    
    
    figure(2);subplot(3,2,2);plot(1:M, correct(:,128),'b-');title('corrected image (:,128)');
    figure(2);subplot(3,2,4);plot(1:M, correct(128,:),'b-');title('corrected image (128,:)');
    figure(2);subplot(3,2,5);plot(1:M, biasSurface(:,128),'b-');title('corrected surface (:,128)');
    figure(2);subplot(3,2,6);plot(1:M, biasSurface(128,:),'b-');title('corrected surface (128,:)');

    figure(3);subplot(2,2,2),mesh(biasSurface);title('Estimated Surface');
    figure(3);subplot(2,2,3);imagesc(img);title('Original Image');
    figure(3);subplot(2,2,4);imagesc(correct);title('Corrected Image');
    
    [label] = E_step(reCorrect',model);
    

end

function biasSurface = biasSurf(img, biasParams)
    [M,N,P] = size(img);
    [X, Y] = meshgrid(1:M, 1:N);
    biasSurface = biasParams(1)*X.^2 + biasParams(2)*Y.^2 + biasParams(3)*X.*Y + biasParams(4)*X + biasParams(5)*Y + biasParams(6);;
end

function model = M_step(correct,label)
correctedImage = reshape(correct,[size(correct,1)*size(correct,2) size(correct,3)]);

nk = sum(label,1);
w = nk/size(label,1);

mu(1) = sum(label(:,1) .* correctedImage) /sum(label(:,1));
mu(2) = sum(label(:,2) .* correctedImage) /sum(label(:,2));

Sigma(:,:,1) = sum(label(:,1) .* (correctedImage - mu(1)).^2) / sum(label(:,1));
Sigma(:,:,2) = sum(label(:,2) .* (correctedImage - mu(2)).^2) / sum(label(:,2));


model.mu = mu;
model.Sigma = Sigma;
model.w = w;
end

function [label] = E_step(X, model)

mu = model.mu;
Sigma = model.Sigma;
w = model.w;

n = size(X,2);
k = 2;
label = zeros(n, k);
norm_label = zeros(1,65536);
for i = 1:k
   norm_label = norm_label + w(i)*(1/sqrt(2*pi*Sigma(:,:,i)))*exp((-(X-mu(i)).^2)/(2*Sigma(:,:,i)));
end

for i = 1:k
label(:,i) = w(i)*(1/sqrt(2*pi*Sigma(:,:,i)))*exp((-(X-mu(i)).^2)/(2*Sigma(:,:,i)));
label(:,i) = label(:,i)./norm_label';
end
end

function biasParams = getC(model, label, reImage,img)
    mu = model.mu;
    sigma = model.Sigma;
    wi1 = label(:,1) / sigma(:,:,1);
    wi2 = label(:,2) / sigma(:,:,2);
    wi = wi1 + wi2;
    W = wi;
    nomt = (wi1 * mu(1) + wi2 * mu(2));
    domt = wi;
    residue = reImage' - (nomt ./ domt);
    [M,N,P] = size(img);
    [X, Y] = meshgrid(1:M, 1:N);
    x = reshape(X, [M*N,1]);
    y = reshape(Y, [M*N,1]);
    A = [x.*x, y.*y, x.*y, x, y, ones(M*N,1)];
    columeA = W.*A;
    columeR = W.*residue;
    biasParams = inv(transpose(A) * columeA) * transpose(A) * columeR;
end
function model=init_values(IStruct)

param = IStruct.param;
mu(1) = param(1);
mu(2) = param(4);
w(1) = param(3);
w(2) = param(6);
sigma(:,:,1) = param(2);
sigma(:,:,2) = param(5);

model.mu = mu;
model.Sigma = sigma;
model.w = w;
end



