% -------------------------------------------------------------------------
% �Գ����ӵĵ������ȷֲ�
% M,N: x,y�������ߵ�Ԫ�ĸ���;
% lambda: ��Ų�����
% dx,dy��x,y����Ԫ�ļ��
% ImΪ��һ���������ֵ
% -------------------------------------------------------------------------

function Ideal_Pattern(Amp,M,N,lambda,dx,dy,Im,theta)%��������ƽ�淽��ͼ
midM=(M-1)/2;midN=(N-1)/2;
k=2*pi/lambda;%����
%% E�淽��ͼ
num = numel(theta);
fEtheta=zeros(1,num);

for m=1:1:M
    for n=1:1:N
        fEtheta=fEtheta+Amp(m,n)*exp(1i*k*(m-midM-1)*dx.*sin(theta));
    end
end
fEtheta=60*Im*abs(cos(0.5*pi*sin(theta))./cos(theta) ).*abs(fEtheta);
figure(2);
plot(180*theta/pi,20*log10(fEtheta./max(fEtheta)));%ֱ�������ͼ�ҹ�һ������
hold on;
title('��һ��E�����۷���ͼ');xlabel('theta');ylabel('f_Etheta');

figure(3);
% polarplot(theta,20*log(max(fEtheta)./fEtheta));%�������һ������
polarplot(theta,20*log(fEtheta));
pax=gca;
pax.ThetaDir='clockwise';
pax.ThetaZeroLocation='top';
title('E�淽��ͼ���۽��');

%% H�淽��ͼ
fHtheta=0;
for m=1:M
    for n=1:N
        fHtheta=fHtheta+Amp(m,n)*exp(1i*k*(n-midN-1)*dy.*sin(theta));
    end
end
fHtheta=60*Im*abs(fHtheta);
figure(4);
plot(180*theta/pi,20*log10(fHtheta./max(fHtheta)));%ֱ�������ͼ�ҹ�һ������
hold on;
title('��һ��H�����۷���ͼ');xlabel('theta');ylabel('f_Htheta');

figure(5);
% polarplot(theta,20*log(max(fEtheta)./fEtheta));%�������һ������
polarplot(theta,20*log(fHtheta));

pax=gca;
pax.ThetaDir='clockwise';
pax.ThetaZeroLocation='top';
title('H�淽��ͼ���۽��');

disp('__________Complete the calculation of ideal pattern__________');
end