function [Ey] = Ey_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im,Mc,Nc,deltax,deltay)
k=2*pi/lambda;%波长
midM=(M-1)/2;midN=(N-1)/2;
midMc=(Mc-1)/2;midNc=(Nc-1)/2;
Xm=[];Yn=[];%第(m, n)号阵列单元的中心点的坐标
Xm1=[];Xm11=[]; %第(m,n)号阵列单元的上下臂的横坐标
Ey=0;%初始化结果
%% 
for n=1:1:N
    Yn(n)=(n-midN-1)*dy;
end

% 在这里Xm1(m)是阵子上臂的坐标，验证的思路是第22个阵子的上臂坐标应该是0.025
for m=1:1:M
    Xm1(m)=(m-midM-1)*dx+lambda/4;
end
% 在这里Xm11(m)是阵子下臂的坐标，验证的思路是第22个阵子的下臂坐标应该是-0.025
for m=1:1:M
    Xm11(m)=(m-midM-1)*dx-lambda/4;
end

% 计算Yl2
for n=1:1:Nc
    Yl2(n)=(n-midNc-1)*deltay;
end

% midMc=(Mc-1)/2;midNc=(Nc-1)/2;% Mc1和Nc2为Mc和Nc的中位数
% XL1=(L1-midMc-1)*deltax;
% YL2=(L2-midNc-1)*deltay;
%% 

for m=1:1:M
    for n=1:1:N
        R1mn=sqrt((Xm1(m)-XL1)^2+(Yn(n)-YL2)^2+d^2);
        R2mn=sqrt((Xm11(m)-XL1)^2+(Yn(n)-YL2)^2+d^2);
        Rhon=sqrt((YL2-Yn).^2+d^2);
        Ey=Ey+Amp(m,n)*((XL1-Xm1(m))*(exp(-1i*k*R1mn))/(R1mn*Rhon(n))+((XL1-Xm11(m))*exp(-1i*k*R2mn))/(R2mn*Rhon(n))); 
    end
end        
end
        
        
