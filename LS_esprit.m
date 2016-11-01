%LS_ESPRIT ALOGRITHM
%DOA ESTIMATION BY LS_ESPRIT ALOGRITHM
clear all;
%close all;
clc;

source_number=2;%信元数
sensor_number=8;%原阵元数
m=7;%子阵元数
N_x=1024; %信号长度
snapshot_number=N_x;%快拍数
w=[pi/4 pi/6].';%信号频率
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2))/2;%信号波长  
d=0.5*l;%阵元间距

snr=0;%信噪比/dB
source_doa=[45 60];%两个信号的入射角度
A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l)].';%阵列流型

s=sqrt(10.^(snr/10))*exp(j*w*[0:N_x-1]);%仿真信号
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号

x1=x(1:m,:);%子阵1接受的数据矢量
x2=x(2:(m+1),:);%子阵2接受的数据矢量

%对两个子阵的模型进行合并
X=[x1;x2];
R=X*X'/snapshot_number;

%对R进行奇异值分解
[U,S,V]=svd(R);
R=R-S(2*m,2*m)*eye(2*m);
[U,S,V]=svd(R);
Us=U(:,1:source_number);
disp(Us);
Us1=Us(1:m,:);
Us2=Us((m+1):2*m,:);


%按照公式得到旋转不变矩阵M
M=pinv(Us1)*Us2;
disp('M');
disp(M);
%对得到的旋转不变矩阵进行特征分解
[Vm,Dm]=eig(M);
disp(Dm);
Dm=(diag(Dm)).';
doa=-asin(angle(Dm)/pi)*180/pi;
disp(doa);