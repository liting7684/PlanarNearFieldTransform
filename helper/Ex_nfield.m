function [Ex] = Ex_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im)
k=2*pi/lambda;%����
midM=(M-1)/2;midN=(N-1)/2;
Xm=[];Yn=[];%��(m, n)�����е�Ԫ�����ĵ������
Xm1=[];Xm11=[]; %��(m,n)�����е�Ԫ�����±۵ĺ�����
Ex=0;%��ʼ�����
%% 
for n=1:1:N
    Yn(n)=(n-midN-1)*dy;
end

% ������Xm1(m)�������ϱ۵����꣬��֤��˼·�ǵ�22�����ӵ��ϱ�����Ӧ����0.025
for m=1:1:M
    Xm1(m)=(m-midM-1)*dx+lambda/4;
end
% ������Xm11(m)�������±۵����꣬��֤��˼·�ǵ�22�����ӵ��±�����Ӧ����-0.025
for m=1:1:M
    Xm11(m)=(m-midM-1)*dx-lambda/4;
end

% midMc=(Mc-1)/2;midNc=(Nc-1)/2;% Mc1��Nc2ΪMc��Nc����λ��
% XL1=(L1-midMc-1)*deltax;
% YL2=(L2-midNc-1)*deltay;
%% 

for m=1:1:M
    for n=1:1:N
        R1mn=sqrt((Xm1(m)-XL1)^2+(Yn(n)-YL2)^2+d^2);
        R2mn=sqrt((Xm11(m)-XL1)^2+(Yn(n)-YL2)^2+d^2);
        Ex=Ex+Amp(m,n)*(exp(-1i*k*R1mn)/R1mn+exp(-1i*k*R2mn)/R2mn); 
    end
end

Ex=-1j*30*Im*Ex;
        
end
        
        
