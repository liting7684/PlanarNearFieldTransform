function NFtrans_Algorithm(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,d,Mc,Nc)
%% --------------------------------------------------
%1.�����Զ���任�ﵽ��E�淽��ͼ����

L1=Mc*deltax*sin(theta)/lambda;
Ex=0;f_Etheta_NF=0;%��ʼ��
for m=0:Mc-1
    for n=0:Nc-1
        XL1=(m-Mc/2)*deltax;
        YL2=(n-Nc/2)*deltay;
        Ex=Ex_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im);%����һ��Ĳ����糡ֵ
        f_Etheta_NF=f_Etheta_NF+Ex*exp(1j*2*pi*m*L1/Mc) ;%����ʽ���㷽��ͼ����
    end
end

f_Etheta_NF=abs(f_Etheta_NF)*deltax*deltay/lambda;%��ʽ5-62



figure(2);
% plot(180*theta/pi,20*log(f_Enf2ff));%��Զ�任����ͼֱ�������ͼ
plot(180*theta/pi,20*log10(f_Etheta_NF./max(f_Etheta_NF)),'-.r');%��Զ�任E�淽��ͼֱ�������ͼ
xlabel('theta');ylabel('f_Etheta_NF');title('��Զ�任E���һ��ֱ�����귽��ͼ');

figure(7);
polarplot(theta,20*log(f_Etheta_NF));%�������ͼ
title('��Զ���任E�淽��ͼ');
pax=gca;
pax.ThetaDir='clockwise';
pax.ThetaZeroLocation='top';


%% --------------------------------------------------
%2.�����Զ���任�ﵽ��H�淽��ͼ����

L2=Nc*deltay*sin(theta)/lambda;
Ex=0;f_Htheta_NF=0;%��ʼ��
for m=0:Mc-1
    for n=0:Nc-1

        XL1=(m-Mc/2)*deltax;
        YL2=(n-Nc/2)*deltay;
        Ex=Ex_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im);%����һ��Ĳ����糡ֵ
        f_Htheta_NF=f_Htheta_NF+Ex*exp(1j*2*pi*n*L2/Nc) ;%����ʽ���㷽��ͼ����
    end
end

f_Htheta_NF=abs(f_Htheta_NF).*abs(cos(theta))*deltax*deltay/lambda;%��ʽ5-62


figure(4);
% plot(180*theta/pi,20*log(f_Enf2ff));%��Զ�任����ͼֱ�������ͼ
plot(180*theta/pi,20*log10(f_Htheta_NF./max(f_Htheta_NF)),'-.r');%��Զ�任E�淽��ͼֱ�������ͼ
xlabel('theta');ylabel('f_Htheta_NF');title('��Զ�任H���һ��ֱ�����귽��ͼ');

figure(9);
polarplot(theta,20*log(f_Htheta_NF));%�������ͼ
title('��Զ���任H�淽��ͼ');
pax=gca;
pax.ThetaDir='clockwise';
pax.ThetaZeroLocation='top';

disp('__________Complete the calculation of pattern by NFtransition Algorithm__________');
end





