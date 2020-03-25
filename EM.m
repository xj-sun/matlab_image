
% implemation of estimating the parameters (mean, standard deviation, and class
%probabilities) for mixtures of Gaussians with EM

%%
close all; clear;


file=input('Please input your file: ', 's');
IStruct=load(file);
data = IStruct.testima;



[M,N,P]=size(data);
len=M*N;
k=3;
d=1;
X=zeros(1, len);
X(1,:)=data(:);
    
 % display original image, histgram, density curve
subplot(2,2,1);
imshow(data); 
subplot(2,2,2);
imhist(data); 
XX=data(:);
bin=0:255;
[counts, centers]=hist(XX, bin);
f1=counts/sum(counts);
subplot(2,2,3)
plot(f1)     

Method=input('Please choose your method(0=EM, 1=fminsearch):');

if (Method==0)
    
    model= init_values();
 
    max_iter = 250; % iterater steps
   
    n=size(X,2);

    [label] = E_step(X,model);

    for iter = 2:max_iter
   
        range=0:255;
        pa = [model.w(1),model.w(2),model.w(3),model.mu(1),model.mu(2),model.mu(3),model.Sigma(1),model.Sigma(2),model.Sigma(3)];

        
        f2=pa(1)*(1/(sqrt(2*pi*pa(7))))*exp((-(range-pa(4)).^2)/(2*pa(7)))+ pa(2)*(1/(sqrt(2*pi*pa(8))))*exp((-(range-pa(5)).^2)/(2*pa(8)))+pa(3)*(1/(sqrt(2*pi*pa(9))))*exp((-(range-pa(6)).^2)/(2*pa(9)));
        f2=f2/sum(f2);
      
        subplot(2,2,4);
        plot(f2)
        pause(0.07) 
        
        model = M_step(X,label);

        [label] = E_step(X,model);
        
    end
    
elseif (Method == 1)
    model= init_values();
    
    pa = [model.w(1),model.w(2),model.mu(1),model.mu(2),model.mu(3),model.Sigma(1),model.Sigma(2),model.Sigma(3)];
    
    optima=fminsearch(@func_fit, pa,[],f1);
    
    op_pa= [optima(1), optima(2),optima(3),optima(4),optima(5),optima(6),optima(7),optima(8)];

    
    range=0:255;

    f2=op_pa(1)*(1/(sqrt(2*pi*op_pa(6))))*exp((-(range-op_pa(3)).^2)/(2*op_pa(6)))+ op_pa(2)*(1/(sqrt(2*pi*op_pa(7))))*exp((-(range-op_pa(4)).^2)/(2*op_pa(7)))+(1-op_pa(1)- op_pa(2))*(1/(sqrt(2*pi*op_pa(8))))*exp((-(range-op_pa(5)).^2)/(2*op_pa(8)));
    f2=f2/sum(f2);
    
    subplot(2,2,4);
    plot(f2)
    
end

function var = func_fit(pa, histo)

    range=0:255;
    f=pa(1)*(1/(sqrt(2*pi*pa(6))))*exp((-(range-pa(3)).^2)/(2*pa(6)))+ pa(2)*(1/(sqrt(2*pi*pa(7))))*exp((-(range-pa(4)).^2)/(2*pa(7)))+(1-pa(1)- pa(2))*(1/(sqrt(2*pi*pa(8))))*exp((-(range-pa(5)).^2)/(2*pa(8)));
    f=f/sum(f);
    var=sum((f-histo).^2);
   
end

function model=init_values()

    init_mu=input('Please type your means or press 0 to run with default([30,60,170]): ');
    init_sig=input('Please type your means or press 0 to run with default([10,10,10]): ');
    init_w=input('Please type your means or press 0 to run with default([0.3,0.3,0.4]): ');
    
    if init_mu==0
        mu = [30,60,170];
        model.mu=mu;
    else
        mu=init_mu;
        model.mu=mu;
    end
    
    if init_sig==0
        for i = 1:3
            Sigma(:,:,i) = 10;
        end
        model.Sigma = Sigma;
    else
        for i=1:3
            Sigma(:,:,i)=init_sig(i);
        end
        model.Sigma=Sigma;
    end
    
    if init_w==0
        w=[0.3,0.3,0.4];
        model.w = w;  
    else
        w=init_w;
        model.w=w;
    end


end


function [label] = E_step(X, model)

   mu = model.mu;
   Sigma = model.Sigma;
   w = model.w;

   n = size(X,2);
   k = 3;
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

function model = M_step(X, label)

   k = size(label,2);
   nk = sum(label,1);

   w = nk/size(label,1);

   mu = bsxfun(@times, X*label, 1./nk);

   Sigma = zeros(1,1,k);
   delta = 1e-8;
   r= sqrt(label);
   for i = 1:k
       Xo = bsxfun(@minus,X,mu(:,i));
       Xo = bsxfun(@times,Xo,r(:,i)');
       Sigma(:,:,i) = Xo*Xo'/nk(i) + delta;
   end
   rate = Sigma/sum(Sigma);
   for i = 1:k
       if rate(i) > 0.5
          Sigma(:,:,i) = rand(1)*min(Sigma);
       end
   end

   model.Sigma = Sigma;
   model.mu = mu;
   model.w = w;
   
end
