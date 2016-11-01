%RB_ESPRIT ALOGRITHM
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

%构造权矩阵W
%ac=zeros(sensor_number,sensor_number);
%for k=0:(sensor_number-1)
%gama(k+1)=k*2*pi/sensor_number;
%ac(:,k+1)=exp(-j*((sensor_number-1)/2)*(gama(k+1)/pi))*exp(j*(0:sensor_number-1)'*(gama(k+1)/pi));
%end
%W=ac;
W=[exp(-j*(sensor_number-1)/2*2*0*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*0*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*1*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*1*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*2*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*2*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*3*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*3*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*4*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*4*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*5*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*5*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*6*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*6*pi/sensor_number) exp(-j*(sensor_number-1)/2*2*7*pi/sensor_number)*exp(j*(0:sensor_number-1)'*2*7*pi/sensor_number)];
disp(W);
%将阵列接收的数据从复数转换成实数
Y=W'*x;
Y1=[real(Y) imag(Y)];
disp(Y1);
R1=Y1*Y1'/(2*snapshot_number);
[U,S,V]=svd(R1);
disp(U);
Es=U(:,1:source_number);
disp(Es);
KS1=[1 cos(pi/sensor_number) cos(2*pi/sensor_number) cos(3*pi/sensor_number) cos(4*pi/sensor_number) cos(5*pi/sensor_number) cos(6*pi/sensor_number) cos(7*pi/sensor_number)];
KS1=diag(KS1);
for i=1:sensor_number-1
    KS1(i,i+1)=cos(i*pi/sensor_number);
end
KS1(sensor_number,1)=(-1).^sensor_number;
KS2=[0 sin(pi/sensor_number) sin(2*pi/sensor_number) sin(3*pi/sensor_number) sin(4*pi/sensor_number) sin(5*pi/sensor_number) sin(6*pi/sensor_number) sin(7*pi/sensor_number)];
KS2=diag(KS2);
for ii=1:sensor_number-1
    KS2(ii,ii+1)=sin(ii*pi/sensor_number);
end

M=pinv(KS1*Es)*(KS2*Es);
%对得到的旋转不变矩阵进行特征分解
[Vm,Dm]=eig(M);
disp(Dm);
beta(1)=-2*atan(Dm(1,1));
beta(2)=-2*atan(Dm(2,2));
disp(beta);
%反解得到信号到达角
doa(1)=asin(beta(1)/pi)*180/pi;
doa(2)=asin(beta(2)/pi)*180/pi;
disp(doa);