function phase_retrieval(Ex1_amplitude,Ex2_amplitude,Ex_data1,lambda,k,Mc,Nc,midMc,midNc,deltax,deltay,d1,d2)
%--------------------------------------------------------------------------
theta0=-pi+2*pi*rand(Mc,Nc);%产生随机初始相位
% theta0=angle(Ex_data1);
%--------------------------------------------------------------------------
kz=zeros(Mc,Nc);
%求解kz
for L1=-Mc/2:1:Mc/2-1
    for L2=-Nc/2:1:Nc/2-1
        flag=k^2-(2*pi*L1/Mc/deltax)^2-(2*pi*L2/Nc/deltay)^2;%判断kz^2-kx^2-ky^2是否大于零
        if  flag>0
            kz(L1+Mc/2+1,L2+Nc/2+1)=sqrt(k^2-(2*pi*L1/Mc/deltax)^2-(2*pi*L2/Nc/deltay)^2);
        else
            kz(L1+Mc/2+1,L2+Nc/2+1)=-1j*sqrt((2*pi*L1/Mc/deltax)^2+(2*pi*L2/Nc/deltay)^2-k^2);
        end
    end
end

for L1=0:1:Mc-1
    for L2=0:1:Nc-1

        kz(L1+1,L2+1)=sqrt(k^2-(2*pi*L1/Mc/deltax)^2-(2*pi*L2/Nc/deltay)^2);

    end
end
%--------------------------------------------------------------------------

Num_iter=500;
% 这个地方存在问题
Ex1_iter=Ex1_amplitude.*exp(1j*theta0);
for i=1:Num_iter
    temp1=fftshift(ifft2(Ex1_iter,Mc,Nc));
    Ex2_iter=fft2(temp1.*exp(1j*kz*(d1-d2)),Mc,Nc);%用1平面场强求出2平面的场强
    Ex2_angle=angle(Ex2_iter);%重新构造2平面的场强，相位为迭代相位
    Ex2_iter=Ex2_amplitude.*exp(1j*Ex2_angle);%幅值为真实幅值
    temp2=fftshift(ifft2(Ex2_iter,Mc,Nc));
    Ex1_iter=fft2(temp2.*exp(-1j*kz*(d1-d2)),Mc,Nc);%用2平面场强求出1平面的场强
    error=sum(sum((abs(Ex1_iter)-Ex1_amplitude).^2))/sum(sum(Ex1_amplitude.^2));%误差计算
    disp(error);
    Ex1_angle=angle(Ex1_iter);%重新构造1平面的场强，相位为迭代相位
    Ex1_iter=Ex1_amplitude.*exp(1j*Ex1_angle);%幅值为真实幅值
end

%近远场变换程序  直接用公式计算了
theta=linspace(-pi/4,pi/4,700);%theta角度范围精度
L1=Mc*deltax*sin(theta)/lambda;
Ex=0;f_Etheta_NF=0;%初始化
for m=0:Mc-1
    for n=0:Nc-1
        XL1=(m-Mc/2)*deltax;
        YL2=(n-Nc/2)*deltay;
        Ex=Ex1_iter(m+1,n+1);
        f_Etheta_NF=f_Etheta_NF+Ex*exp(1j*2*pi*m*L1/Mc) ;%带公式计算方向图函数

    end
end

f_Etheta_NF=abs(f_Etheta_NF)*deltax*deltay/lambda;%公式5-62
figure(2);
% plot(180*theta/pi,20*log(f_Enf2ff));%近远变换方向图直角坐标绘图
plot(180*theta/pi,20*log10(f_Etheta_NF./max(f_Etheta_NF)),'-.r');%近远变换E面方向图直角坐标绘图
xlabel('theta');ylabel('f_Etheta_NF');title('近远变换E面归一化直角坐标方向图');

end


