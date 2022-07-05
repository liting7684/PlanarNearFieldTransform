% -------------------------------------------------------------------------
% MATLAB Code 
% Generated on: 07-Mar-2022 15:38:42

% Generate:The main function;
% Function: generate the antenna pattern using NearField data and comparing
% with ideal pattern;
% Current_Excitation(M,N,dx,dy);计算电流分布Amn
% Ideal_Pattern(Amp,M,N,lambda,dx,dy,Im,theta);计算理想平面方向图
% NFtrans_Algorithm(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,d,Numc);近远场变换计算方向图

% variables:
% M,N: x,y方向天线单元的个数;
% lambda: 电磁波波长
% Im: 半波振子波腹电流
% dx,dy: x,y方向天线单元的间距
% theta: theta角度范围精度
% deltax,deltay: 在x轴和y轴上的采样间距
% d: 采样平面距离阵列平面的间距
% Mc,Nc: 计算机仿真x,y方向的采样点数目,取2的倍数
% -------------------------------------------------------------------------

clear all;close all;clc;
M=43; N=31;
% freq = 9375e6;
% lambda = physconst('LightSpeed')/freq;
lambda=32;
Im=1;%半波振子波腹电流
dx=0.7*lambda;dy=0.7*lambda;%x,y方向天线单元的间距
theta=linspace(-pi/4,pi/4,600);%theta角度范围精度
phi=linspace(0,2*pi,600);%phi角度范围精度
deltax=0.45*lambda;deltay=0.45*lambda;% deltax和deltay为在x轴和y轴上的采样间距
d=4*lambda;%采样平面距离阵列平面的间距
Mc=128; Nc=128;%计算机仿真x,y方向的采样点数目,取2的倍数
padding = 1;%充0个数
disp('__________this is the beginning__________')

Amp=Current_Excitation(M,N,dx,dy);%计算电流分布Amn
%%
for i = 11:1:15
    for j = 18:1:21
        Amp(i,j)=0;
    end
end

surf(Amp)
%%
Ideal_Pattern(Amp,M,N,lambda,dx,dy,Im,theta)%计算理想平面方向图

disp('__________plotting the pattern by NFtransition Algorithm__________');
NFtrans_Algorithm(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,d,Mc,Nc);%近远场变换计算方向图
disp('__________plotting the pattern by NFtransition Algorithm Using FFT__________');
data_nf2ff = nf2ff_planar_fft(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,phi,d,Mc,Nc,padding);%用平面波谱推导远场方向图
disp('__________plotting phi=0 and phi=90 cut__________');
plotFFPhiCut(data_nf2ff,[0,pi/2]);
disp('__________this is the end__________')

