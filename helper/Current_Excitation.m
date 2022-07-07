function [Amp]=Current_Excitation(M,N,dx,dy)%计算电流分布Amn

midM=(M-1)/2;midN=(N-1)/2;
Xm=[];Yn=[];%第(m, n)号阵列单元的中心点的坐标
for m=1:1:M
    Xm(m)=(m-midM-1)*dx;
end
% 在这里yn(n)是阵子中心y轴的坐标，验证的思路是第16个阵子为0
for n=1:1:N
    Yn(n)=(n-midN-1)*dy;
end
ampx=ones(1,M);%令x方向为等幅度
ampy=ones(1,N);
%y方向余弦分布
for i=1:N
    ampy(i)=cos(pi*(i-midN-1)/N);
end

% x方向切比雪夫分布
SLL=-55;%副瓣电平
ampx=Chebyshev_Dist_Odd(M,SLL);
% ampx=chebwin1(43,-55);

Amp=ampx'*ampy;
[ymesh,xmesh]=meshgrid(Yn,Xm);
figure(1);
surf(xmesh,ymesh,Amp);
xlabel('xdimension');ylabel('ydimension');title('对称阵子幅度分布');
disp('__________Complete the calculation of current__________');
% disp(Amp);
end