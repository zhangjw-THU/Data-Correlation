close all
clear

%% ����5֮

%% ����1
f = 2000000;%��������Ƶ��
a = 1*10^12;
tau = 5*10^(-6);

Time = 2*tau;
fs = 50*10^6;%����Ƶ��
N = fs*Time;  %���������
dt = 1/fs;    %����ʱ����
T = (0:N-1)*dt;%ÿ���������ʱ��ڵ�

x1 = sin(2*pi*f*T).*exp((-a*(T-tau).^2)/2);

figure(1)
plot(T,x1,'b')
xlabel('t/s');
title('������˳���źŲ���');

%% ����2
delta = [205*10^(-9),210*10^(-9),215*10^(-9)];
x2 = sin(2*pi*f*(T-delta(1))).*exp((-a*(T-delta(1)-tau).^2)/2);
x3 = sin(2*pi*f*(T-delta(2))).* exp((-a*(T-delta(2)-tau).^2)/2);
x4 = sin(2*pi*f*(T-delta(3))).*exp((-a*(T-delta(3)-tau).^2)/2);
figure(2)
plot(T,x1,'b');
hold on
plot(T,x2,'y');
xlabel('t/s');
title('������˳��-�����źŲ���');
legend('˳���ź�','������ʱ205ns')

figure(3)
plot(T,x1,'b');
hold on
plot(T,x3,'c');
xlabel('t/s');
title('������˳��-�����źŲ���');
legend('˳���ź�','������ʱ210ns')

figure(4)
plot(T,x1,'b');
hold on
plot(T,x4,'r');
xlabel('t/s');
title('������˳��-�����źŲ���');
legend('˳���ź�','������ʱ215ns');

figure(5)
plot(T,x1,'b');
hold on
plot(T,x2,'y');
hold on
plot(T,x3,'c');
hold on
plot(T,x4,'r');
xlabel('t/s');
title('������˳��-�����źŲ���');
legend('˳���ź�','������ʱ205ns','������ʱ210ns','������ʱ215ns');


%% ����3
% �����ݽ��в�ֵ
fs_0 = 1*10^9;
dt_0 = 1/fs_0;
N_0 = fs_0*Time;
T_0 = (0:N_0-1)*dt_0;
x1_0 = x1;
x2_0 = x2;
x1 = interp1(T,x1,T_0,'pchip');
x2 = interp1(T,x2,T_0,'pchip');
x3 = interp1(T,x3,T_0,'pchip');
x4 = interp1(T,x4,T_0,'pchip');

figure(20)
plot(T_0,x1,'b');
hold on
plot(T_0,x2,'y');
xlabel('t/s');
title('������˳��-�����źŲ���');
legend('˳���ź�','������ʱ205ns')


[a,b] = xcorr(x1,x2);
figure(6)
plot(b*(1/fs_0),a,'y');
xlabel('t/s');
title('�ź�����Է���');
[Rmax,i] = max(a);
N1 = (i-length(x1))*(1/fs);

[a,b] = xcorr(x1,x3);
figure(7)
plot(b*(1/fs_0),a,'c');
xlabel('t/s');
title('�ź�����Է���');
[Rmax,i] = max(a);
N2 = (i-length(x1))*(1/fs);

[a,b] = xcorr(x1,x4);
figure(8)
plot(b*(1/fs_0),a,'r');
xlabel('t/s');
title('�ź�����Է���');
[Rmax,i] = max(a);
N3 = (i-length(x1))*(1/fs);


%% ����4������
NoiseVar(1) = 0;
for i= 1:10000
    
    Noise1 = NoiseVar(i)*rand(1,length(x1_0));
    Noise2 = NoiseVar(i)*rand(1,length(x2_0));
    
    x1_0_noised = x1_0 + Noise1;
    x2_0_noised = x2_0 + Noise2;
    
    x1_0_noised_itp = interp1(T,x1_0_noised,T_0,'pchip');
    x2_0_noised_itp = interp1(T,x2_0_noised,T_0,'pchip');
    
    [a,b] = xcorr(x1_0_noised_itp,x2_0_noised_itp);
    [Rmax,I] = max(a);
    Res(i) = (I-length(x1_0_noised_itp))*(1/fs_0);
    NoiseVar(i+1) =  NoiseVar(i)+0.0001;
end
figure(9)
plot( NoiseVar(1:i),Res)
title('�����Խ����Ӱ��');
xlabel('����ǿ��');
ylabel('���(t/s)');



