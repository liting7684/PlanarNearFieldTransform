%% 
%--------------------------------------------------------------------------
% MATLAB Code
% Generated on: 29-June-2022 13:17:42
% Author: THYLOVEZJ 
% Function: Use two planar amplitude scanner plane data to retrieve phase,and use
% nf2ff technology to push near-field to far-field
% variables:
% M,N: x,y方向天线单元的个数;
% lambda: 电磁波波长
% Im: 半波振子波腹电流
% dx,dy: x,y方向天线单元的间距
% theta: theta角度范围精度
% deltax,deltay: 在x轴和y轴上的采样间距
% d: 采样平面距离阵列平面的间距
% Mc,Nc: 计算机仿真x,y方向的采样点数目,取2的倍数
%--------------------------------------------------------------------------
%% 
disp('__________Set up first scanner plane parameter__________')
close all;clc;
% load first amplitude scanner plane data 
% The space between first scanner plane and aperture is 4*lambda
M=43; N=31;
% freq = 9375e6;
% lambda = physconst('LightSpeed')/freq;
lambda=32;
Im=1;%半波振子波腹电流
dx=0.7*lambda;dy=0.7*lambda;%x,y方向天线单元的间距
theta=linspace(-pi/4,pi/4,600);%theta角度范围精度
phi=linspace(0,pi/2,600);%phi角度范围精度
deltax=0.45*lambda;deltay=0.45*lambda;% deltax和deltay为在x轴和y轴上的采样间距
d1=4*lambda;%采样平面距离阵列平面的间距
Mc=128; Nc=128;%计算机仿真x,y方向的采样点数目,取2的倍数
midMc=Mc/2;midNc=Nc/2;
k=2*pi/lambda;
%%
disp('__________Calculate first plane two orthogonal electric field.__________')
Amp=Current_Excitation(M,N,dx,dy);%计算第一个面的电流分布Amn
%计算第一个面的Ex和Ey
L1=Mc*deltax*sin(theta)/lambda;
Ex1=0;f_Etheta_NF=0;%初始化
for m=0:Mc-1
    for n=0:Nc-1
        XL1=(m-Mc/2)*deltax;
        YL2=(n-Nc/2)*deltay;
        Ex1(m+1,n+1)=Ex_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d1,Im);%计算Ex的采样电场值
    end
end
for m=0:Mc-1
    for n=0:Nc-1
        XL1=(m-Mc/2)*deltax;
        YL2=(n-Nc/2)*deltay;
        Ey1(m+1,n+1)=Ey_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d1,Im,Mc,Nc,deltax,deltay);%计算Ey的采样电场值
    end
end
%%
disp('__________Set up second scanner plane parameter__________')
d2=4.8*lambda;
for m=0:Mc-1
    for n=0:Nc-1
        XL1=(m-Mc/2)*deltax;
        YL2=(n-Nc/2)*deltay;
        Ex2(m+1,n+1)=Ex_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d2,Im);%计算Ex的采样电场值
    end
end

Ey2=0;
for m=0:Mc-1
    for n=0:Nc-1
        XL1=(m-Mc/2)*deltax;
        YL2=(n-Nc/2)*deltay;
        Ey2(m+1,n+1)=Ey_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d2,Im,Mc,Nc,deltax,deltay);%计算Ey的采样电场值
    end
end
disp('__________Calculate phase retrieval Algorithm__________')
phase_retrieval(abs(Ex1),abs(Ex2),angle(Ex1),lambda,k,Mc,Nc,midMc,midNc,deltax,deltay,d1,d2)
disp('__________This is the end__________')
