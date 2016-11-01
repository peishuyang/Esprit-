%LS_ESPRIT ALOGRITHM
%DOA ESTIMATION BY LS_ESPRIT ALOGRITHM
clear all;
%close all;
clc;

bbb=zeros(1,11);
for kk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%信噪比(dB)

aaa=zeros(1,300);
for k=1:300

source_number=1;%信元数
sensor_number=8;%原阵元数
m=7;%子阵元数
N_x=1024; %信号长度
snapshot_number=N_x;%快拍数
w=pi/4;%信号频率
l=2*pi*3e8/w;%信号波长 
d=0.5*l;%阵元间距

source_doa=50;%信号的入射角度

A=[exp(-j*(0:sensor_number-1)*d*2*pi*sin(source_doa*pi/180)/l)].';
s=10.^((snr(kk)/2)/10)*exp(j*w*[0:N_x-1]);%仿真信号

%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号

x1=x(1:m,:);%子阵1接受的数据矢量
x2=x(2:(m+1),:);%子阵2接受的数据矢量

%对两个子阵的模型进行合并
X=[x1;x2];
R=X*X'/snapshot_number;
%对R进行奇异值分解
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
estimated_doa=-asin(angle(Dm(1,1))/pi)*180/pi;
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%求解均方根误差和标准偏差
%E_doa=sum(aaa(1,:))/300;%做300次试验的均值
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%做300次试验的均方根误差
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kk)=RMSE_doa;
end
disp(bbb);
plot(snr,bbb(1,:),'k*-');
save LS_ESPRIT_snr_rmse.mat;
%TLS_ESPRIT
bbb=zeros(1,11);
for kk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%信噪比(dB)

aaa=zeros(1,300);
for k=1:300
    
s=sqrt(10.^(snr(kk)/10))*exp(j*w*[0:N_x-1]);%仿真信号
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号

x1=x(1:m,:);%子阵1接受的数据矢量
x2=x(2:(m+1),:);%子阵2接受的数据矢量

%对两个子阵的模型进行合并
X=[x1;x2];
R=X*X'/snapshot_number;
%对R进行奇异值分解
[U,S,V]=svd(R);
Us=U(:,1:source_number);
disp(Us);
Us1=Us(1:m,:);
Us2=Us((m+1):2*m,:);
Us12=[Us1 Us2];
%形成矩阵Us12
Us12=[Us1,Us2];
%对“Us12'*Us12”进行特征分解，得到矩阵E
[E,Sa,Va]=svd(Us12'*Us12);
disp('E');
disp(E);
disp(Sa);
%将E分解为四个小矩阵
E11=E(1,1);
E12=E(1,2);
E21=E(2,1);
E22=E(2,2);
%按照公式得到旋转不变矩阵M
M=-(E12*(inv(E22)));
disp('M');
disp(M);
%对得到的旋转不变矩阵进行特征分解
[Vm,Dm]=eig(M);
disp(Dm);
estimated_doa=-asin(angle(Dm(1,1))/pi)*180/pi;
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%求解均方根误差和标准偏差
%E_doa=sum(aaa(1,:))/300;%做300次试验的均值
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%做300次试验的均方根误差
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kk)=RMSE_doa;
end
disp(bbb);
hold on
plot(snr,bbb(1,:),'rs-');
save TLS_ESPRIT_snr_rmse.mat;
%TAM
bbb=zeros(1,11);
for kk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%信噪比(dB)

aaa=zeros(1,300);
for k=1:300
    
s=sqrt(10.^(snr(kk)/10))*exp(j*w*[0:N_x-1]);%仿真信号
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号

R=x*x'/snapshot_number;
disp(R);
%对接收到的数据协方差矩阵进行奇异值分解,分解得到信号子空间及大特征值
[U,S,V]=svd(R);
Us=U(:,1:source_number);
Ss=S(1:source_number,1:source_number);
disp(Ss);
Vs=V(:,1:source_number);
%利用旋转不变子空间思想构造矩阵B
B=Us*(Ss^(1/2));
B1=B(1:(sensor_number-1),:);
B2=B(2:sensor_number,:);
%提取Us的子向量构造Us1和Us2
%Us1=Us(1:(sensor_number-1),:);
%Us2=Us(2:sensor_number,:);
%利用以上关系得到最小二乘解
%D=pinv(Us1*((Ss)^(1/2)))*Us2*((Ss)^(1/2));
D=pinv(B1)*B2;
%对D进行特征分解，由特征值可得到对应的N个信号的到达角
[Vd,Dd]=eig(D);
disp(Dd);
estimated_doa=-asin(angle(Dd(1,1))/pi)*180/pi;
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%求解均方根误差和标准偏差
%E_doa=sum(aaa(1,:))/300;%做300次试验的均值
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%做300次试验的均方根误差
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kk)=RMSE_doa;
end
disp(bbb);
hold on
plot(snr,bbb(1,:),'bd-');
save TAM_snr_rmse.mat
%R_ESPRIT
bbb=zeros(1,11);
for kkk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%信噪比(dB)

aaa=zeros(1,300);
for k=1:300
    
s=sqrt(10.^(snr(kkk)/10))*exp(j*w*[0:N_x-1]);%仿真信号

%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+j*randn(sensor_number,N_x));%加了高斯白噪声后的阵列接收信号
%构造矩阵Qm和Q2L
J=[0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0];%四阶置换矩阵
Qm=(1/sqrt(2))*[eye(4) j*eye(4);J -j*J];
JJ=zeros(snapshot_number,snapshot_number);
for ii=1:snapshot_number
      JJ(ii,snapshot_number+1-ii)=1;
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
%反解得到信号到达角
estimated_doa =real(asin(beta(1)/pi)*180/pi); %计算来波方向
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%求解均方根误差和标准偏差
%E_doa=sum(aaa(1,:))/300;%做300次试验的均值
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%做300次试验的均方根误差
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kkk)=RMSE_doa;
end
disp(bbb);
hold on
plot(snr,bbb(1,:),'gx-');
save resprit_snr_rmse.mat
%RB_ESPRIT
bbb=zeros(1,11);
for kkk=1:11
snr=[-10 -8 -6 -4 -2 0 2 4 6 8 10];%信噪比(dB)

aaa=zeros(1,300);
for k=1:300
    
s=sqrt(10.^(snr(kkk)/10))*exp(j*w*[0:N_x-1]);%仿真信号

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
%反解得到信号到达角
estimated_doa =real(asin(beta(1)/pi)*180/pi); %计算来波方向
disp(estimated_doa);

aaa(:,k)=estimated_doa;

disp('estimated_doa');
disp(estimated_doa);
end
disp(aaa);

%求解均方根误差和标准偏差
%E_doa=sum(aaa(1,:))/300;%做300次试验的均值
E_doa=mean(aaa);
disp(E_doa);

RMSE_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%做300次试验的均方根误差
disp('RMSE_doa');
disp(RMSE_doa);

bbb(:,kkk)=RMSE_doa;
end
disp(bbb);
hold on
plot(snr,bbb(1,:),'mo-');
save RB_esprit_snr_rmse.mat

legend('LS-ESPRIT','TLS-ESPRIT','TAM','实值空间ESPRIT','实值波束空间ESPRIT');
xlabel('信噪比（snr）/dB');
ylabel('估计均方根误差');
grid on;


