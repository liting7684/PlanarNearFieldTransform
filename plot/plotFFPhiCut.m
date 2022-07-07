function plotFFPhiCut(data_ff,phi)
%--------------------------------------------------------------------------
% 找到phi等于所求phi=0切面的索引
i = find(data_ff.phi==phi(1,1));
% 通过该索引找到phiCut等于固定值时theta的坐标
theta=data_ff.theta(i);
% 找到该索引对应的Etheta的最大值
maxValue = max(max(data_ff.Etheta(i)));
% nf2ff_cut即是phi等于所给值时Etheta
nf2ff_cut =  abs(data_ff.Etheta(i))/maxValue;
%--------------------------------------------------------------------------
figure(14)
title("phi=0时Etheta和Ephi")
plot(theta*180/pi,20*log10(nf2ff_cut),'r--');
hold on
%--------------------------------------------------------------------------
% 找到该索引对应的Ephi的最大值
% maxValue = max(max(data_ff.Ephi(i)));
% nf2ff_cut即是phi等于所给值时Ephi
nf2ff_cut =  abs(data_ff.Ephi(i))/maxValue;
%--------------------------------------------------------------------------
title("phi=0时Ephi")
plot(theta*180/pi,20*log10(nf2ff_cut),'m-.');
legend("phi=0时Etheta","phi=0时Ephi")
%--------------------------------------------------------------------------
% nf2ff_cut即是phi等于所给值时Eabs
nf2ff_cut =  abs(data_ff.Eabs(i))/maxValue;
%--------------------------------------------------------------------------
figure(15)
title("phi=0时Eabs")
plot(theta*180/pi,20*log10(nf2ff_cut),'r--');
%--------------------------------------------------------------------------
% 找到phi等于所求phi=90切面的索引

i = find(data_ff.phi==phi(1,2));
% 通过该索引找到phiCut等于固定值时theta的坐标
theta=data_ff.theta(i);
% 找到该索引对应的Etheta的最大值
maxValue = max(max(data_ff.Ephi(i)));
% nf2ff_cut即是phi等于所给值时Etheta
nf2ff_cut =  abs(data_ff.Etheta(i))/maxValue;
%--------------------------------------------------------------------------
figure(16)
title("phi=pi/2时Etheta和Ephi")
plot(theta*180/pi,20*log10(nf2ff_cut),'r--');
hold on
%--------------------------------------------------------------------------
% nf2ff_cut即是phi等于所给值时Ephi
nf2ff_cut =  abs(data_ff.Ephi(i))/maxValue;
%--------------------------------------------------------------------------
title("phi=pi时Ephi")
plot(theta*180/pi,20*log10(nf2ff_cut),'m-.');
legend("phi=pi/2时Etheta","phi=pi/2时Ephi")
%--------------------------------------------------------------------------
% nf2ff_cut即是phi等于所给值时Eabs
nf2ff_cut =  abs(data_ff.Eabs(i))/maxValue;
%--------------------------------------------------------------------------
figure(17)
title("phi=pi/2时Eabs")
plot(theta*180/pi,20*log10(nf2ff_cut),'r--');
%--------------------------------------------------------------------------
end