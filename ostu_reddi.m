%implementation of the 2-class Otsu thresholding algorithm and of the 3-class Reddi
%thresholding algorithm

file=input('What is your file name?', 's')
IStruct=load(file)
I=IStruct.testima

Method=input('which method do you like (0=Ostu, 1=Reddi)')




%%


% Ostu
if (Method==0)
    subplot(2,2,1)
    imshow(I); % show the original one
    subplot(2,2,2)
    
    hold on
    
    imhist(I); % show hisg

    n=imhist(I); 
    N=sum(n); 
    
    for i=1:256
        P(i)=n(i)/N; % value of probability
    end
    
    max=0; % start point(maximaize the value)
    
    for T=2:255      
        p0=sum(P(1:T)); 
        p1=sum(P(T+1:256)); 
        u0=dot([0:T-1],P(1:T))/p0;
        u1=dot([T:255],P(T+1:256))/p1; 
        k=p0*p1*((u1-u0)^2); % the value of curent threshold

        if k>max 
            max=k; 
            threshold=T-1; % update threshold
        end
    end
    
    plot([threshold, threshold], [0,2000],'red')
    
    hold off
    
    res=im2bw(I,threshold/255); % Convert the original image by our threshold
    subplot(2,2,3)
    imshow(res); % show the final result

% Reddi
elseif  (Method==1)
    subplot(2,2,1)
    imshow(I); % show the  Original one
    subplot(2,2,2)
    
    hold on
    
    imhist(I); % show the hisg
    n=imhist(I); % Compute the histogram
    N=sum(n); % sum the values of all the histogram values
    max=0; % start point(maximaize the value)

    for i=1:256
        P(i)=n(i)/N; % value of probability
    end


    th1=86 %  256 / 3
    th2=170 % 2*th1

    x1 = th1;
    y1 = 0:255;
    x2 = th2;
    y2 = 0:255;

    plot([x1,x1], [0,2000], 'red')
    plot([x2,x2], [0,2000], 'red')
    
    p0=sum(P(1:th1)); 
    p1=sum(P(th1+1:th2)); 
    p2=sum(P(th2+1:256)); 
    u0=dot([0:th1-1],P(1:th1))/p0; % mean of class0
    u1=dot([th1:th2-1],P(th1+1:th2))/p1; % mean of class1
    u2=dot([th2:255],P(th2+1:256))/p2; % mean of class2
    los1=(u0+u1)/2-th1;
    los2=(u1+u2)/2-th2;
    count = 0
    % start point
    
    while count < 20 && (abs(los1)>1 || abs(los2)>1)
        th1=th1+round(los1);
        th2=th2+round(los2);

        p0=sum(P(1:th1)); 
        p1=sum(P(th1+1:th2)); 
        p2=sum(P(th2+1:256)); 
        u0=dot([0:th1-1],P(1:th1))/p0; 
        u1=dot([th1:th2-1],P(th1+1:th2))/p1; 
        u2=dot([th2:255],P(th2+1:256))/p2; 
        los1=(u0+u1)/2-th1;
        los2=(u1+u2)/2-th2;

        x1 = th1;
        x2 = th2;
        plot([x1,x1], [0,2000], 'red')
        plot([x2,x2], [0,2000], 'red')
        pause(0.3)
        count=count + 1

    end
    
    hold off
   
    th=[th1;th2]
    res = imquantize(I,th); % convert image by our threshold
    subplot(2,2,4)
    imshow(res,[]); % show the result
end
