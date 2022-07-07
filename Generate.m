% -------------------------------------------------------------------------
% MATLAB Code 
% Generated on: 07-Mar-2022 15:38:42

% Generate:The main function;
% Function: generate the antenna pattern using NearField data and comparing
% with ideal pattern;
% Current_Excitation(M,N,dx,dy);��������ֲ�Amn
% Ideal_Pattern(Amp,M,N,lambda,dx,dy,Im,theta);��������ƽ�淽��ͼ
% NFtrans_Algorithm(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,d,Numc);��Զ���任���㷽��ͼ

% variables:
% M,N: x,y�������ߵ�Ԫ�ĸ���;
% lambda: ��Ų�����
% Im: �벨���Ӳ�������
% dx,dy: x,y�������ߵ�Ԫ�ļ��
% theta: theta�Ƕȷ�Χ����
% deltax,deltay: ��x���y���ϵĲ������
% d: ����ƽ���������ƽ��ļ��
% Mc,Nc: ���������x,y����Ĳ�������Ŀ,ȡ2�ı���
% -------------------------------------------------------------------------

clear all;close all;clc;
M=43; N=31;
% freq = 9375e6;
% lambda = physconst('LightSpeed')/freq;
lambda=32;
Im=1;%�벨���Ӳ�������
dx=0.7*lambda;dy=0.7*lambda;%x,y�������ߵ�Ԫ�ļ��
theta=linspace(-pi/4,pi/4,600);%theta�Ƕȷ�Χ����
phi=linspace(0,2*pi,600);%phi�Ƕȷ�Χ����
deltax=0.45*lambda;deltay=0.45*lambda;% deltax��deltayΪ��x���y���ϵĲ������
d=4*lambda;%����ƽ���������ƽ��ļ��
Mc=128; Nc=128;%���������x,y����Ĳ�������Ŀ,ȡ2�ı���
padding = 1;%��0����
disp('__________this is the beginning__________')

Amp=Current_Excitation(M,N,dx,dy);%��������ֲ�Amn
%%
for i = 11:1:15
    for j = 18:1:21
        Amp(i,j)=0;
    end
end

surf(Amp)
%%
Ideal_Pattern(Amp,M,N,lambda,dx,dy,Im,theta)%��������ƽ�淽��ͼ

disp('__________plotting the pattern by NFtransition Algorithm__________');
NFtrans_Algorithm(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,d,Mc,Nc);%��Զ���任���㷽��ͼ
disp('__________plotting the pattern by NFtransition Algorithm Using FFT__________');
data_nf2ff = nf2ff_planar_fft(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,phi,d,Mc,Nc,padding);%��ƽ�沨���Ƶ�Զ������ͼ
disp('__________plotting phi=0 and phi=90 cut__________');
plotFFPhiCut(data_nf2ff,[0,pi/2]);
disp('__________this is the end__________')

