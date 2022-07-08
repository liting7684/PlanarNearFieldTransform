function [] = plotFFPolarization(data_nf2ff,phi_cut)

%plotFFPolarization(data_nf2ff)����Ҫ��theta_range��phi_range�ĵ���һ��
normalized = true;
logarithmic = true;

theta_cut=max(data_nf2ff.theta);
% Find values in cutting plane
i = find(data_nf2ff.phi==phi_cut);
j = find(data_nf2ff.theta==theta_cut);

% Flip one have of the values to make continous plot
Etheta_cut_angles = data_nf2ff.theta(i);
Ephi_cut_angles = data_nf2ff.phi(j);

% %���ڲ���ƽ�治�������Σ����Ϳ�һ��������ʱ��������
% %���öԲ������ٵı߽��в�ֵ�ķ�����ʹ�ó��Ϳ�һ��
% Etheta_cut_angles = interp2(x,y,Etheta_cut',xq,yq,'spline');

Etheta_maxValue = max(max(data_nf2ff.Etheta(i)));
Ephi_maxValue = max(max(data_nf2ff.Ephi(j)));

if normalized == true
    %ȡ��ȡabsͼ�񼸺�û�в��
    %Etheta_cut =  data_nf2ff.Etheta(i)/Etheta_maxValue;
	%Ephi_cut =  data_nf2ff.Ephi(j)/Ephi_maxValue;
    Etheta_cut =  abs(data_nf2ff.Etheta(i))/Etheta_maxValue;
    Ephi_cut =  abs(data_nf2ff.Ephi(j))/Ephi_maxValue;
else
    %Etheta_cut =  data_nf2ff.Etheta(i);
    %Ephi_cut =  data_nf2ff.Ephi(j);
    Etheta_cut =  abs(data_nf2ff.Etheta(i));
    Ephi_cut =  abs(data_nf2ff.Ephi(j));
end


Eref=sin(Ephi_cut_angles).*Etheta_cut+cos(Ephi_cut_angles).*Ephi_cut;
Ecross=cos(Ephi_cut_angles').*Etheta_cut-sin(Ephi_cut_angles').*Ephi_cut;
Ecross=Ecross(:,1);



if logarithmic == true
    figure
    title("this is Eref and Ecross!");
    plot(Etheta_cut_angles*180/pi,20*log10(Eref),Etheta_cut_angles*180/pi,20*log10(Ecross));
    legend("Ecross","Eref");
    xlabel('theta');
    ylabel('E_{ref}��E_{cross}');
%     xlabel('theta');
%     ylabel('E_{cross}');
elseif logarithmic == false
    figure
    plot(Etheta_cut_angles*180/pi,Eref,Etheta_cut_angles*180/pi,Ecross);
    legend("Ecross","Eref");
    xlabel('theta');
    ylabel('E_{ref}��E_{cross}');
%     xlabel('theta');
%     ylabel('E_{cross}');
end
hold on
end
