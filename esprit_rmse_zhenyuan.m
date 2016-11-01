%LS_ESPRIT ALOGRITHM
%DOA ESTIMATION BY LS_ESPRIT
clear all;
%close all;
clc;

bbb=zeros(1,10);
for kk=1:10
sensor_number=[3 4 6 8 10 12 14 16 18 20];%阵元数
m=[2 3 5 7 9 11 13 15 17 19];%子阵元数
aaa=zeros(1,300);
for k=1:300

source_number=1;%信元数
N_x=1024; %信号长度
snapshot_number=N_x;%快拍数
w=pi/4;%信号频率
l=2*pi*3e8/w;%信号波长  
d=0.5*l;%阵元间距
snr=-2;%信噪比(dB)
source_doa=50;%信号的入射角度
A=[exp(-j*(0:sensor_number(kk)-1)*d*2*pi*sin(source_doa*pi/180)/l)].';

s=10.^(snr/20)*exp(j*w*[0:N_x-1]);%仿真信号
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number(kk),N_x)+j*randn(sensor_number(kk),N_x));%加了高斯白噪声后的阵列接收信号

x1=x(1:m(kk),:);%子阵1接受的数据矢量
x2=x(2:(m(kk)+1),:);%子阵2接受的数据矢量

%对两个子阵的模型进行合并
X=[x1;x2];
R=X*X'/snapshot_number;
%对R进行奇异值分解
[U,S,V]=svd(R);
Us=U(:,1:source_number);
disp(Us);
Us1=Us(1:m(kk),:);
Us2=Us((m(kk)+1):2*m(kk),:);

%按照公式得到旋转不变矩阵M
M=pinv(Us1)*Us2;
disp('M');
disp(M);
%对得到的旋转不变矩阵进行特征分解
[Vm,Dm]=eig(M);
disp(Dm);
estimated_source_doa=-asin(angle(Dm(1,1))/pi)*180/pi;
disp(estimated_source_doa);

aaa(:,k)=estimated_source_doa;
end
disp(aaa);

%求解均方根误差和标准偏差
%E_source_doa=sum(aaa(1,:))/300;%做300次试验的均值
E_source_doa=mean(aaa);
disp(E_source_doa);

RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%做300次试验的均方根误差
disp('RMSE_source_doa');
disp(RMSE_source_doa);

bbb(:,kk)=RMSE_source_doa;
end
disp(bbb);

plot(sensor_number,bbb(1,:),'k*-');
save LS_ESPRIT_zhenyuan_RMSE.mat;
%TLS_ESPRIT
bbb=zeros(1,10);
for kk=1:10
sensor_number=[3 4 6 8 10 12 14 16 18 20];%阵元数
m=[2 3 5 7 9 11 13 15 17 19];%子阵元数

aaa=zeros(1,300);
for k=1:300
A=[exp(-j*(0:sensor_number(kk)-1)*d*2*pi*sin(source_doa*pi/180)/l)].';

s=10.^(snr/20)*exp(j*w*[0:N_x-1]);%仿真信号
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number(kk),N_x)+j*randn(sensor_number(kk),N_x));%加了高斯白噪声后的阵列接收信号
x1=x(1:m(kk),:);%子阵1接受的数据矢量
x2=x(2:(m(kk)+1),:);%子阵2接受的数据矢量

%对两个子阵的模型进行合并
X=[x1;x2];
R=X*X'/snapshot_number;
%对R进行奇异值分解
[U,S,V]=svd(R);
Us=U(:,1:source_number);
disp(Us);
Us1=Us(1:m(kk),:);
Us2=Us((m(kk)+1):2*m(kk),:);
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
estimated_source_doa=-asin(angle(Dm(1,1))/pi)*180/pi;

aaa(:,k)=estimated_source_doa;
end
disp(aaa);

%求解均方根误差和标准偏差
E_source_doa=sum(aaa(1,:))/300;%做300次试验的均值
disp('E_source_doa');

RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%做300次试验的均方根误差
disp('RMSE_source_doa');
disp(RMSE_source_doa);

bbb(:,kk)=RMSE_source_doa;
end
disp(bbb);
hold on
plot(sensor_number,bbb(1,:),'rs-');
save TLS_ESPRIT_zhenyuan_rmse.mat;
%TAM ALOGRITHM
bbb=zeros(1,10);
for kk=1:10
sensor_number=[3 4 6 8 10 12 14 16 18 20];%阵元数

aaa=zeros(1,300);
for k=1:300
A=[exp(-j*(0:sensor_number(kk)-1)*d*2*pi*sin(source_doa*pi/180)/l)].';
s=10.^(snr/20)*exp(j*w*[0:N_x-1]);%仿真信号
%x=awgn(s,snr);
x=A*s+(1/sqrt(2))*(randn(sensor_number(kk),N_x)+j*randn(sensor_number(kk),N_x));%加了高斯白噪声后的阵列接收信号

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
B1=B(1:(sensor_number(kk)-1),:);
B2=B(2:sensor_number(kk),:);
%提取Us的子向量构造Us1和Us2
%Us1=Us(1:(sensor_number(kk)-1),:);
%Us2=Us(2:sensor_number(kk),:);
%利用以上关系得到最小二乘解
%D=pinv(Us1*((Ss)^(1/2)))*Us2*((Ss)^(1/2));
D=pinv(B1)*B2;
%对D进行特征分解，由特征值可得到对应的N个信号的到达角
[Vd,Dd]=eig(D);
disp(Dd);
estimated_source_doa=-asin(angle(Dd(1,1))/pi)*180/pi;
%计算信号到达方向角
aaa(:,k)=estimated_source_doa;
end
disp(aaa);

%求解均方根误差和标准偏差
E_source_doa=sum(aaa(1,:))/300;%做300次试验的均值
disp('E_source_doa');

RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%做300次试验的均方根误差
disp('RMSE_source_doa');
disp(RMSE_source_doa);

bbb(:,kk)=RMSE_source_doa;
end
disp(bbb);

hold on
plot(sensor_number,bbb(1,:),'bd-');

save TAM_zhenyuan_RMSE.mat;

legend('LS-ESPRIT','TLS-ESPRIT','TAM');
xlabel('阵元数目');
ylabel('估计均方根误差');
grid on;





