function NFtrans_Algorithm(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,d,Mc,Nc)
%% --------------------------------------------------
%1.计算近远场变换达到的E面方向图函数

L1=Mc*deltax*sin(theta)/lambda;
Ex=0;f_Etheta_NF=0;%初始化
for m=0:Mc-1
    for n=0:Nc-1
        XL1=(m-Mc/2)*deltax;
        YL2=(n-Nc/2)*deltay;
        Ex=Ex_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im);%计算一点的采样电场值
        f_Etheta_NF=f_Etheta_NF+Ex*exp(1j*2*pi*m*L1/Mc) ;%带公式计算方向图函数
    end
end

f_Etheta_NF=abs(f_Etheta_NF)*deltax*deltay/lambda;%公式5-62



figure(2);
% plot(180*theta/pi,20*log(f_Enf2ff));%近远变换方向图直角坐标绘图
plot(180*theta/pi,20*log10(f_Etheta_NF./max(f_Etheta_NF)),'-.r');%近远变换E面方向图直角坐标绘图
xlabel('theta');ylabel('f_Etheta_NF');title('近远变换E面归一化直角坐标方向图');

figure(7);
polarplot(theta,20*log(f_Etheta_NF));%极坐标绘图
title('近远场变换E面方向图');
pax=gca;
pax.ThetaDir='clockwise';
pax.ThetaZeroLocation='top';


%% --------------------------------------------------
%2.计算近远场变换达到的H面方向图函数

L2=Nc*deltay*sin(theta)/lambda;
Ex=0;f_Htheta_NF=0;%初始化
for m=0:Mc-1
    for n=0:Nc-1

        XL1=(m-Mc/2)*deltax;
        YL2=(n-Nc/2)*deltay;
        Ex=Ex_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im);%计算一点的采样电场值
        f_Htheta_NF=f_Htheta_NF+Ex*exp(1j*2*pi*n*L2/Nc) ;%带公式计算方向图函数
    end
end

f_Htheta_NF=abs(f_Htheta_NF).*abs(cos(theta))*deltax*deltay/lambda;%公式5-62


figure(4);
% plot(180*theta/pi,20*log(f_Enf2ff));%近远变换方向图直角坐标绘图
plot(180*theta/pi,20*log10(f_Htheta_NF./max(f_Htheta_NF)),'-.r');%近远变换E面方向图直角坐标绘图
xlabel('theta');ylabel('f_Htheta_NF');title('近远变换H面归一化直角坐标方向图');

figure(9);
polarplot(theta,20*log(f_Htheta_NF));%极坐标绘图
title('近远场变换H面方向图');
pax=gca;
pax.ThetaDir='clockwise';
pax.ThetaZeroLocation='top';

disp('__________Complete the calculation of pattern by NFtransition Algorithm__________');
end





