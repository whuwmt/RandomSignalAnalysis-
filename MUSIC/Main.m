% ====== 一维超分辨传统MUSIC算法 ====== %
% =======......2021.06.10......====== %  
clear all;clc;close all;

%%%%%%% ==== ..........参数初始化
f0 = 7.65e6;           %载波频率
w0 = 2*pi*f0;          %载波角频率
c = 3e8;               %光速
Lambda = c/f0;         %载频波长
Snap = 1000;            %快拍数
Snr = 1000;               %信噪比
% ==== 天线阵列位置
AtnX = [-46,0,-18,-36,-54,-72,-49.9099998474121,-58.0299987792969];
AtnY = [29.1700000762939,0,0,0,0,0,14.5200004577637,-14.3599996566772];
scatter(AtnX,AtnY);
[~,NumAtn] = size(AtnX);%阵列数

% ==== 数学模型构建:非相干信号源
Doa  = [30 50 70];     %要估计的角度
NumSource = length(Doa);%信源数
A = A_theta(AtnX,AtnY,NumAtn,w0,Doa);
Nr = zeros(NumAtn,Snap); %加性噪声
% === 非相干信号
Signal = randn(NumSource,Snap)+1j*randn(NumSource,Snap);
X= Snr*A*Signal+Nr;
% ===== 非相干信号源(LFM)
fs = 1e3;
u1 = 1;
u2 = 1.5;
u3 = 2;
t = (0:Snap-1)/fs;
s1=exp(j*2*pi*(f0*t+1/2*u1*t.^2));%产生的线性调频信号1包络(点序列)                                     
s2=exp(j*2*pi*(f0*t+1/2*u2*t.^2));%产生的线性调频信号2包络(点序列)  
s3=exp(j*2*pi*(f0*t+1/2*u3*t.^2));%产生的线性调频信号1包络(点序列)
S = [s1;s2;s3];
X = A*S;
% ===== MUSIC估计
flag =1;        %采用特征值分解
Research = [0 90 1];%开始搜索角度，结束搜索角度，步长
[P,theta] = MUSIC_Tranditon(X,NumSource,NumAtn,Snap,AtnX,AtnY,flag,w0,Research);
% ====== plot,画图
plot(theta*180/pi,10*log10(P/max(P)));
axis([-180 180 -50 +10])
xlabel('入射角/（度）','fontsize',10)
ylabel('Pmusic','fontsize',10)
title('MUSIC非相干信源一维估计')