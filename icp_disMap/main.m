% implementation two rigid-body registration methods. One based on a distance map
%and the other on the Iterative Closest Point algorithm.

%prompt = 'Please input the mat:\n';
 %path = input(prompt,'s');
 %load(path);
 %load path;
load './testhw5_5.mat';

image = zeros(256,256);

if (method == 1)
    PointOptimization(image,xc,yc,pointx,pointy);
end

if(method == 2)
    costFunction = ICP(image,xc,yc,pointx,pointy);
end

function cos = ICP(image,xc,yc,pointx,pointy)


P = [pointx;pointy];

plot(P(:,1),P(:,2),'r*-');
count = 0;
cos = [];
figure(1),title('Registration Process');
imagesc(image);
hold on;

S = {};
for i = 1 : numel(xc) - 1
    S{i,1} = [xc(i),yc(i)];
    S{i,2} = [xc(i+1),yc(i+1)];
    S{i,3} = (xc(i+1) - xc(i))^2 + (yc(i+1) - yc(i))^2;
end

while(1)
    
    matchedPoints = findMatch(S,P);
       
    % break situation and show the registration result
    count = count + 1;
    cos(count) = calDistances(P,matchedPoints);
    if (count > 1)
        
        if (abs(calDistances(P,matchedPoints) - cos(count-1)) < 1)
            figure(2),title('Registration Result');
            imagesc(image);
            hold on;
            plot(xc,yc,'b*-');
            plot(Tpoints(1,:),Tpoints(2,:),'r*-');
            break;
        end
    end
    
    [R, shift] = register(P,matchedPoints);
    %center = [128; 128];
    Tpoints = R * P + repmat(shift,[1,size(P,2)]);
    plot(xc,yc,'b*-');
    plot(Tpoints(1,:),Tpoints(2,:),'r*-');
    pause(0.2);
    P = Tpoints;
       
end

end

function match = findMatch(lineSegments,P)

match= [];
for i = 1:size(P,2)
    x3 = P(1,i);
    y3 = P(2,i);
    nearestPoint = findNearest(lineSegments,x3,y3);
    match(:,i) = nearestPoint;
end

end

function near = findNearest(S,x3,y3)

%intinal loss value
loss = 7777777;

for i = 1:size(S,1)
    x1 = S{i,1}(1);
    y1 = S{i,1}(2);
    x2 = S{i,2}(1);
    y2 = S{i,2}(2);
    
    denomintor = S{i,3};
    
    miu = ((x3 - x1) * (x2 -  x1) + (y2 - y1) * (y3 - y1)) / denomintor;
    if(miu > 0 && miu < 1)
        tempPoint(1) = x1 + miu * (x2 - x1);
        tempPoint(2) = y1 + miu * (y2 - y1);
    end
    
    if(miu < 0 || miu > 1)
        if (calDistance(x3,y3,x1,y1) > calDistance(x3,y3,x2,y2))
            tempPoint = [x2,y2];
        else
            tempPoint = [x1,y1];
        end
    end
    
    dist = calDistance(x3,y3,tempPoint(1),tempPoint(2));
    
    if (dist < loss)
        loss = dist;
        near = tempPoint;
    end
end
end

function dist = calDistance(pointx,pointy,boundx,boundy)
dist = (pointx - boundx)^2+(pointy - boundy)^2;
end

function [R, t] = register(X,Y)

Xc= X -repmat(128,1,size(X,2));
Yc= Y -repmat(128,1,size(X,2));
[U, ~, V] = svd(Xc*Yc');        %singular value decomposition of matrix A
R = V*U';
t = mean(Y,2) -R*mean(X,2);
end

function distance = calDistances(X,Y)
distance = sum(sum((X - Y).^2));
end

                
                

