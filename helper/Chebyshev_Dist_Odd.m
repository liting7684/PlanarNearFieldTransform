function [Ix1]=Chebyshev_Dist_Odd(M,SLL)
% 返回切比雪夫分布归一化电流
% M=43;
midM=(M-1)/2;
% SLL=-55;
R0=10^(SLL/(-20));
a0=cosh(acosh(R0)/(M-1));  
a0=cosh(log(R0+sqrt(R0^2-1))/(M-1));

n=midM;

Ix=zeros(1,n+1);

for k=0:n
    for q=k:1:n
        
        Ix(k+1)=Ix(k+1)+(-1)^(n-q)*a0^(2*q)*(factorial(q+n-1)*(2*n)...
        /(factorial(q-k)*factorial(q+k)*factorial(n-q)));
        
    end
end

% even case;
% for k=1:n
%     for q=k:1:n
%         
%         Ix(k)=Ix(k)+(-1)^(n-q)*a0^(2*q-1)*(factorial(q+n-2)*(2*n-1)/...
% (factorial(q-k)*factorial(q+k-1)*factorial(n-q)));
%         
%     end
% end

Ix=Ix./max(Ix);
Ix1=zeros(1,M);
for i=1:midM
    Ix1(i)=Ix(midM-i+2);
end
for i=midM+1:M
    Ix1(i)=Ix(i-midM);
end

% disp(Ix1');
end