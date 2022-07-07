%%
%--------------------------------------------------------------------------
% Function: This .m file solve the real measurement nf2ff process.
% variables:
% Ex_retrieval: 采样恢复出的Ex
% Ey_retrieval: 采样恢复出的Ey
% Mc: x方向上的采样点数
% Nc: y方向上的采样点数
% 
% author: THYLOVEZJ
% 
%--------------------------------------------------------------------------
%%
% 载入实测数据恢复Ex
data = importdata('./data/4GHz-p1-v.dat');
% 恢复出原来矩阵测试 51*35
Mc = 51;
Nc = 35;
data_reformat = reshape(data(:,1),[Mc,Nc]);
E_aplitude = exp(data_reformat./20);
data_reformat = reshape(data(:,2),[Mc,Nc]);
% 恢复出原来测试的相位 51*35
E_phase = reshape(data_reformat,[Mc,Nc]);
% 恢复出原来测试的电场
Ex_retrieval = E_aplitude.*cosd(E_phase)+1i*E_aplitude.*sind(E_phase);
%%
% 载入实测数据恢复Ey
data = importdata('./data/4GHz-p1-h.dat');
data_reformat = reshape(data(:,1),[Mc,Nc]);
E_aplitude = exp(data_reformat./20);
data_reformat = reshape(data(:,2),[Mc,Nc]);
E_phase = reshape(data_reformat,[Mc,Nc]);
Ey_retrieval = E_aplitude.*cosd(E_phase)+1i*E_aplitude.*sind(E_phase);
%%
% other settings
% we use 4G horn,so lambda was set at 75mm
lambda = 75;
% 采样间隔
deltax = 
%%
nf2ff_planar_real(Ex,Ey,lambda,)
