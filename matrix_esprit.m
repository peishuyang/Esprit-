%MATRIX ESPRIT ALOGRITHM
%DOA ESTIMATION BY MATRIX ESPRIT ALOGRITHM
clear all;
close all;
clc;

source_number=2;%信元数
sensor_number=8;%原阵元数
m=7;%子阵元数
N_x=1024; %信号长度
snapshot_number=N_x;%快拍数
w=[pi/4 pi/6].';%信号频率
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2))/2;%信号波长  
d=0.5*l;%阵元间距
array_distance=d;%两组阵元轴线方向的间距
snr=0;%信噪比

source_doa=[45 60];%两个信号的入射角度

A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l)].';%阵列流型

s=sqrt(10.^(snr/10))*exp(j*w*[0:N_x-1]);%仿真信号
%构建矩阵Z
Z=zeros(m,m);
for i=1:(m-1)
    Z(i+1,i)=1;
end

%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号

x1=x(1:m,:);%子阵1接收的数据
x2=x(2:sensor_number,:);%子阵2接收的数据

R11=x1*x1'/snapshot_number;%子阵1的自相关矩阵
R12=x1*x2'/snapshot_number;%子阵1、2的互相关矩阵

D1=eig(R11);
u2=mean(D1(1:m-source_number)); %估计噪声的方差

%去噪得到子阵1、2的自和互协方差矩阵
C11=R11-u2*eye(m);
C12=R12-u2*Z;

%对C11、C12进行广义特征值分解，得到N个大特征值
D=eig(C11,C12);
D_big=zeros(1,source_number);
[Y L]=sort(abs(abs(D)-1));
D=D(L);
disp(D);
D_big=[D(2) D(1)];
%反解得到信号到达角
doa =-asin((angle(D_big)/pi))*180/pi %计算来波方向






