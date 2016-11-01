%R_ESPRIT ALOGRITHM
%DOA ESTIMATION BY R_ESPRIT
clear all;
close all;
clc;

source_number=2;%信元数
sensor_number=8;%阵元数
N_x=1024; %信号长度
snapshot_number=N_x;%快拍数
w=[pi/4 pi/6].';%信号频率
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2))/2;%信号波长  
d=0.5*l;%阵元间距
snr=0;%信噪比

source_doa=[45 60];%两个信号的入射角度
A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(1)*pi/180)/l);exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa(2)*pi/180)/l)].';%阵列流型

s=sqrt(10.^(snr/10))*exp(j*w*[0:N_x-1]);%仿真信号
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号
%构造矩阵Qm和Q2L
J=[0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0];%四阶置换矩阵
Qm=(1/sqrt(2))*[eye(4) j*eye(4);J -j*J];
JJ=zeros(snapshot_number,snapshot_number);
for i=1:snapshot_number
      JJ(i,snapshot_number+1-i)=1;
end
Q2L=(1/sqrt(2))*[eye(snapshot_number) j*eye(snapshot_number);JJ -j*JJ];
%构造Z矩阵
Jm=[0 0 0 0 0 0 0 1;0 0 0 0 0 0 1 0;0 0 0 0 0 1 0 0;0 0 0 0 1 0 0 0;0 0 0 1 0 0 0 0;0 0 1 0 0 0 0 0;0 1 0 0 0 0 0 0;1 0 0 0 0 0 0 0];
xx=(x').';%对数据矢量取共轭
Z=[x Jm*xx*JJ];
%构造变换矩阵Tx
Tx=Qm'*Z*Q2L;
%通过变换矩阵Tx将阵列接收数据转换到实值空间Rt
Rt=Tx*Tx'/(2*snapshot_number);
%对Rt进行奇异值分解，得到其信号子空间Es
[Ur,Sr,Vr]=svd(Rt);
Es=Ur(:,1:source_number);
%利用LS求解旋转不变性
kk=[0 0 0 0 0 0 0].';
K2=[kk eye(sensor_number-1)];%构造矩阵K2
%构造矩阵Qmm
JJJ=[0 0 1;0 1 0;1 0 0];
Qmm=(1/sqrt(2))*[eye(3) [0 0 0].' j*eye(3);0 0 0 sqrt(2) 0 0 0;JJJ [0 0 0].' -j*JJJ];
%运算得到H1和H2
H1=2*real(Qmm'*K2*Qm);
H2=2*imag(Qmm'*K2*Qm);
%得旋转不变矩阵M
M=pinv(H1*Es)*(H2*Es);

%对得到的旋转不变矩阵进行特征分解
[Vm,Dm]=eig(M);
disp(Dm);
beta(1)=-2*atan(Dm(1,1));
beta(2)=-2*atan(Dm(2,2));
disp(beta);
%反解得到信号到达角
doa(1)=real(asin(beta(1)/pi)*180/pi);
doa(2)=real(asin(beta(2)/pi)*180/pi);
disp(doa);
