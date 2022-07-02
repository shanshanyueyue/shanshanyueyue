%利用FFT提取相位，每帧只取一个chirp，求每个帧的相位
%数据 采样率512  帧长40ms  帧数640  每帧chirp数 1
%FFT点数应该两倍于ADC采样点数

clear all;close all;clc;
%% 雷达参数（使用mmWave Studio默认参数）
c=3.0e8;  
B=4e9;           %调频带宽
K=29982e6;       %调频斜率
T=B/K;           %采样时间
Tc=60e-6;        %chirp总周期
fs=1e4;          %采样率
f0=77e9;         %初始频率
lambda=c/f0;     %雷达信号波长
d=lambda/2;      %天线阵列间距0
n_samples=512;   %采样点数/脉冲
N=1024;          %距离向FFT点数
n_chirps=1;      %每帧脉冲数
n_RX=1;          %RX天线通道数
frame = 640;     %帧数

retVal = readDCA1000("test11.bin");
%retVal = readDCA1000("test8_Raw_0.bin");
cdata = zeros(n_RX,n_chirps*n_samples);
% 生成一个n_RX行，n_chirps*n_samples列的零矩阵
allPhases = (n_chirps);
Index = (n_chirps);
for xx=1:frame
    for row = 1:1
      cdata(row,:) = retVal(row,(xx-1)*n_samples*n_chirps*64+1:((xx-1)*64+1)*n_samples*n_chirps);  % retVal的数据格式是512个采样点连续64个为一个frame
    end    % cdata里面保存了每一个frame里面首个chirp里面的512个采样点数据
    data_radar_1 = reshape(cdata(1,:),n_samples,n_chirps);   % RX1 里面的数据变成n_samples(512)*n_chirps(1)的矩阵
    data_radar=[];            
    data_radar(:,:,1)=data_radar_1;     %三维雷达回波数据g
    range_win = hamming(n_samples);   %加海明窗
    range_profile = [];
    for k=1:n_RX
       for m=1:n_chirps
          %temp=data_radar(:,m,k).*range_win;    %加窗函数
          temp_fft=fft(data_radar(:,m,k),N);    %对每个chirp做N点FFT
          range_profile(:,m,k)=temp_fft;
       end
    end   
end

%%%%%%%%%% 一阶fft原图  %%%%%%%%%%
x=1:1024;
y=abs(range_profile(x,1,1));
figure(1);
plot(x,y);
set(gca,'XTick',0:11:1024);
set(gca,'XTickLabel',(0:11:1024)/22);
xlabel('距离(米)');
ylabel('幅值');
title('幅值-距离图');

%%%%%%%%%% 投票结果图  %%%%%%%%%%
for xx=1:frame
    cdata(1,:) = retVal(1,(xx-1)*n_samples*n_chirps*64+1:((xx-1)*64+1)*n_samples*n_chirps);
    % cdata里面保存了每一个frame里面首个chirp里面的512个采样点数据
    temp_fft2=fft(cdata(1,:),N);  %对每个chirp做N点FFT
    range_profile2(1,:)=temp_fft2;
    range_abs2 = (abs(range_profile2));
    anglediffer(xx,:)=angle(temp_fft2);
end
for i=1:1024
    unwrapPhase1(:,i) = unwrap(anglediffer(:,i));
    unwrapPhase1(:,i)=unwrapPhase1(:,i)-mean(unwrapPhase1(:,i));
end
distance=zeros(1,40);
for i=1:40
    imf=emd(unwrapPhase1(:,i));
    imf1(:,i)=imf(1,:);
    avgPhase(1,i)=mean(abs(imf1(:,i)));
    for j=1:640
        if abs(imf1(j,i))>3.5*avgPhase(1,i)
            distance(1,i)=distance(1,i)+1;
        end
    end
end
figure(2);
plot(1:40,distance);
set(gca,'XTick',0:11:40);
set(gca,'XTickLabel',(0:11:40)/22);
title('最终投票结果');
xlabel('距离(米)');
ylabel('票数');
[n,t]=max(distance);

%%%%%%%%%% 相位-时间图  %%%%%%%%%%
for xx=1:frame
    cdata(row,:) = retVal(row,(xx-1)*n_samples*n_chirps*64+1:((xx-1)*64+1)*n_samples*n_chirps); 
    data_radar_1 = reshape(cdata(1,:),n_samples,n_chirps);   
    data_radar=[];            
    data_radar(:,:,1)=data_radar_1;    
    range_profile = [];
    temp_fft=fft(data_radar(:,1,1),N);  
    range_profile(:,1,1)=temp_fft; 
    range_abs = (abs(range_profile));
    [maxVal,maxIndex] = max(range_abs);
    Index(xx) = maxIndex;
    allPhases(t,xx) = angle(range_profile(t));
end
unwrapPhase = unwrap(allPhases(t,:));
figure(3);
plot(unwrapPhase(1:640));
title('相位-时间图');
set(gca,'XTick',0:10:680);
set(gca,'XTickLabel',(0:10:680)/25);
xlabel('时间(秒)');
ylabel('相位');

%%%%%%%%%% 幅值-频率图  %%%%%%%%%%
dropMinLimit=3; %表示在25.6秒的时间内最少滴的滴数
dropMaxLimit=60;%表示在25.6秒的时间内最多滴的滴数，用于过滤
y = abs(fft(unwrapPhase1(:,t),N));  % 对上式进行 N 点 FFT 计算 ，并取模值
f = (1:1024)*frame/N;    % 转换为频率区间
figure(4);
plot(f(1:N/2),y(1:N/2));
title('频谱图');
xlabel('频率(Hz)');
ylabel('幅度');
maxf=y(1);
mindrop=int8((dropMinLimit)/frame*N);
maxdrop=int8((dropMaxLimit)/frame*N);
[m,index]=max(y(mindrop+1:maxdrop,1));

disp("滴数:"+int8((index+mindrop)/1.6));