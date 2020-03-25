function  [theta,shift] = PointOptimization(image,xcontour,ycontour,pointx,pointy)


BW = roipoly(image,xcontour,ycontour);
edges = edge(BW);
edges = double(edges);
imagesc(edges);

global points
points = [pointx;pointy];
global distanceMap;
distanceMap = bwdist(edges);

figure(1),imagesc(distanceMap);
title('Original Shape')
hold on;
plot(pointx,pointy,'r*-');
hold off;

shift = [0,0];
distance = 7777;
global RT;

%find the local minum from dfferet original angle and save the local minum
for i = 1 : 4
    
    theta = (i - 1) * 90;       
    x0 = [theta shift] /1000;    
    [x,fun] = fminsearch(@calDistance,x0);      
    if (fun < distance)
        distance = fun;
        bestParam = x;
        histRT = RT;
    end
    
    RT = [];
end

figure(2),title('Registration Process');
imagesc(distanceMap);
hold on;

for i = 1: 10:size(histRT,1)  
    
theta = histRT(i,1);
shift = histRT(i,2:3);    
phi = theta;
R = [cos(phi),-sin(phi);sin(phi),cos(phi)];
centerPoint = [128;128];
Tpoints = R * (points - repmat(centerPoint,[1,size(points,2)])) + repmat(centerPoint,[1,size(points,2)])+ repmat(shift',[1,size(points,2)]);
plot(Tpoints(1,:),Tpoints(2,:),'g*-');
pause(0.7);
    
end

figure(3)
imagesc(distanceMap);
hold on;

theta = bestParam(1) * 1000;
shift = bestParam(2:3) * 1000;    
phi = theta;
R = [cos(phi),-sin(phi);sin(phi),cos(phi)];
centerPoint = [128;128];
%centerPoint = mean(points,2);

Tpoints = R * (points - repmat(centerPoint,[1,size(points,2)])) + repmat(centerPoint,[1,size(points,2)])+ repmat(shift',[1,size(points,2)]);

plot(Tpoints(1,:),Tpoints(2,:),'g*-');
plot(xcontour,ycontour,'r*-');
title('Registration Result');

end

function dis = calDistance(param)
global points;
global distanceMap;
global RT;

theta = param(1) * 1000;
shift = param(2:3) * 1000;

RT = [RT;[theta,shift]];
phi = theta;
R = [cos(phi),-sin(phi);sin(phi),cos(phi)];
centerPoint = [128;128];
Tpoints = R * (points - repmat(centerPoint,[1,size(points,2)])) + repmat(centerPoint,[1,size(points,2)])+ repmat(shift',[1,size(points,2)]);
dis = interp2(distanceMap,Tpoints(1,:),Tpoints(2,:));
dis = sum(dis(:));

end


