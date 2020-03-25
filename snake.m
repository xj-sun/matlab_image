close all; clear;

%file = input('What is the name of your file?', 's')
%IStruct = load(file) 
pwd = '/Users/xujuansun/Downloads'
Istruct = load([pwd,'/testhw4_3'])
forceType = Istruct.forcetype
std = Istruct.std
sup = Istruct.support
niter = Istruct.Niter
nsam = Istruct.nsample
bin = Istruct.binim
a = Istruct.alpha
b = Istruct.beta
g = Istruct.gamma
eCoef = Istruct.extcoef
bCoef = Istruct.balcoef
iterGVF = Istruct.itergvf
initContourx = Istruct.initcontourx
initContoury = Istruct.initcontoury

snake(bin, forceType, std, sup, niter, nsam, a, b, g, eCoef, bCoef, iterGVF,initContourx,initContoury)






function [] = snake(bin, forceType, std, sup, niter, nsam, a, b, g, eCoef, bCoef, iterGVF,X,Y)
    
    
    invIA = Finternal_Matrix(a, b, g, nsam) % internal force matrix
    
    [FextX, FextY] = Fexternal(bin, forceType, sup, std, iterGVF) % external force
     
    figure(1); imagesc(bin); hold on


    [X, Y] = resample_contour(X', Y', nsam)
    
    % show the curve evolve process
    for j = 1 : niter
        
        [fbalX, fbalY] = Fballoon(X, Y)   
        fextX = interp2(FextX, X, Y, 'linear')
        fextY = interp2(FextY, X, Y, 'linear')
        X = invIA * (X' + eCoef.*g.*fextX' + bCoef*fbalX')
        Y = invIA * (Y' + eCoef.*g.*fextY' + bCoef*fbalY')
    
        X= [X',X(1)];
        Y =[Y',Y(1)];
        [X, Y] = resample_contour(X, Y, nsam)
        
        for i = 1: numel(X) - 1
            
            line([X(i), X(i+1)], [Y(i), Y(i+1)], 'Color', 'red', 'LineWidth',2)
            plot(X(i), Y(i), 'g*')
            
        end
    %     plot([X; X(1)], [Y; Y(1)], 'g-');  
        pause(0.1)
    end
    hold off
    %   plot final curve
    figure; imshow(bin, [0 0 1; 1 0 0]); hold on
    plot(X, Y, 'g*'); hold off
end

function [FbalX, FbalY] = Fballoon(X, Y)

% clockwise direction curve, Fballoon points to outside of the contour
X = circshift(X, -1, 1) - circshift(X, 1, 1); % X(i+1)-X(i-1)
Y = circshift(Y, -1, 1) - circshift(Y, 1, 1); % Y(i+1)-Y(i-1)
mag = sqrt(X.^2 + Y.^2)
FbalX = Y ./ (mag+1e-10)
FbalY = -X ./ (mag+1e-10)

end

function [FextX, FextY] = Fexternal(im, forceType, sup, std, iterGVF)

if forceType == 1                                             % gradient
    im = double(im);
    im = (im - min(im(:))) ./ (max(im(:)) - min(im(:)));
    im = im - 0.005;
    G = fspecial('gaussian', sup*2+1, std);
    im = imfilter(im, G, 'replicate');
    [FextX, FextY] = gradient(im); 
    
elseif forceType == 2                                         % distance
   
    Dis = bwdist(im)
    [FextX, FextY] = gradient(-Dis)
    
elseif forceType == 3                                          % GVF
    
    [FextX, FextY] = GVF(im, iterGVF) 
end
 
figure, quiver(FextX, FextY); axis('image', 'on', 'ij')
end

function invIA = Finternal_Matrix(a, b, g, nsam)

a1 = b
a2 = -a - 4*b
a3 = 2*a + 6*b
a4 = -a - 4*b
a5 = b
a = [a3, a4, a5, zeros(1,nsam-5), a1, a2]
A = zeros(nsam, nsam)
for i = 1:nsam
    A(i,:) = circshift(a, i-1, 2);
end
invIA = inv(eye(nsam) + g.*A);
end

function [u, v] = GVF(im, iterGVF)

% normalize f to the range [0,1]
mu = 0.1
[r, c] = size(im);
im = padarray(im, [1, 1], 'replicate', 'both')
[fx, fy] = gradient(im)   
mag = fx.^2 + fy.^2

% initialize GVF to the gradient
u = fx
v = fy
for i = 1 : iterGVF
  u = u + mu.*4.*del2(u) - (u - fx).*mag
  v = v + mu.*4.*del2(v) - (v - fy).*mag
end

u = u(2:(r+1), 2:(c+1))
v = v(2:(r+1), 2:(c+1))
end

function [x, y] = resample_contour(X, Y, nsam)

x = X
y=Y
len = numel(x)
l=zeros(1, len)
for i =1 : len-1
    l(i + 1) = l(i) + sqrt((x(i+1)-x(i)).^2 + (y(i+1)-y(i)).^2)
end
delta =l(end)/nsam

interpx = interp1(l, x, (0:delta:l(end)))
interpy = interp1(l, y, (0:delta:l(end)))
x = interpx(1:nsam)
y = interpy(1:nsam)
end