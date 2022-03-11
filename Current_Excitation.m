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
for i=1:N
    ampy(i)=cos(pi*(i-midN-1)/N);%y方向余弦分布
end
SLL=-20;
R0=10^(SLL/(-20));
a0=cosh(acosh(R0)/(M-1));  
 
Ix=[];
for j=1:1:midM
    
    Ix(1)=(-1)^(midM-j-1)*a0^(2*(j-1))*(factorial(1+midM-2)*(2*midM)...
        /(2*factorial(j-1)*factorial(j+1-2)*factorial(midM-j+1)));

end

for i=2:midM
    for j=1:1:midM
        
        Ix(i)=(-1)^(midM-j-1)*a0^(2*(j-1))*(factorial(i+midM-2)*(2*midM)...
        /(1*factorial(j-i)*factorial(j+i-2)*factorial(midM-j+1)));
        
    end
end


Amp=ampx'*ampy;
[ymesh,xmesh]=meshgrid(Yn,Xm);
figure(1);
surf(xmesh,ymesh,Amp);
xlabel('xdimension');ylabel('ydimension');title('对称阵子幅度分布');
disp('__________Complete the calculation of current__________');
end