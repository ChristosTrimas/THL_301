clear all;
close all;
clc;

%a.1
%like exercise 1 and 2
N = 500;
b = (sign(randn(4*N,1))+1)/2;

%a.2
%A=1 to get the same bits_to_4_PAM with prvious exercises
A = 1;
X = bits_to_4_PAM(b,A);

%a.3
Xi = X(1:N);
Xq = X(N+1:2*N);

%a.4
T = 0.01;
over = 10;
Ts = T/over;
Fs = 1/Ts;

%we use some values from previous exercises to simulate SRRC
A1 = 4;
a = 0.5;
N = 2048;

%Frequency
F = -Fs/2:Fs/N:Fs/2-Fs/N;

[phi,t] = srrc_pulse(T,Ts,A1,a);

%signals
Xi_delta = Fs * upsample(Xi,over);
Xi_delta_conv = conv(Xi_delta,phi)*Ts;

Xq_delta = Fs * upsample(Xq,over);
Xq_delta_conv = conv(Xq_delta,phi)*Ts;

%time
t1 = 0:Ts:((N+N*(over-1))-1)*Ts;
Ti_conv = linspace(t(1)+t1(1), t(end) + t1(end),length(Xi_delta_conv)); 
Tq_conv = linspace(t(1)+t1(1), t(end) + t1(end),length(Xq_delta_conv));

figure;
subplot(2,1,1);
plot(Ti_conv,Xi_delta_conv);
title('Xi(t) signal');
xlabel('time');
ylabel('Xi(t)');

subplot(2,1,2);
plot(Tq_conv,Xq_delta_conv);
title('Xq(t) signal');
xlabel('time');
ylabel('Xq(t)');

%fasmatiki piknothta isxios
PxXi = (abs(fftshift(fft(Xi_delta_conv,N))).^2)./length(Ti_conv)*Ts;
PxXq = (abs(fftshift(fft(Xq_delta_conv,N))).^2)./length(Tq_conv)*Ts;

figure;
subplot(2,1,1)
plot(F,PxXq);
title('Periodogram of PxXq');
xlabel('F(Hz)');
ylabel('PxXq');

subplot(2,1,2);
plot(F,PxXi);
title('Periodogram of PxXi');
xlabel('F(Hz)');
ylabel('PxXi');

figure;
subplot(2,1,1)
semilogy(F,PxXq);
title('Periodogram of PxXq with semilogy');
xlabel('F(Hz)');
ylabel('PxXq');

subplot(2,1,2);
semilogy(F,PxXi);
title('Periodogram of PxXi with semilogy');
xlabel('F(Hz)');
ylabel('PxXi');

%a.5
F0 = 200;
Xi_mod = 2*Xi_delta_conv.*cos(2*pi*F0*Ti_conv);
Xq_mod = -2*Xq_delta_conv.*sin(2*pi*F0*Tq_conv);

figure;
subplot(2,1,1);
plot(Tq_conv,Xq_mod);
title('Xq_mod');
xlabel('Xq_mod');
ylabel('time');

subplot(2,1,2);
plot(Ti_conv,Xi_mod);
title('Xi_mod');
xlabel('Xi_mod');
ylabel('time');


PxXq_mod = (abs(fftshift(fft(Xq_mod,N))).^2)./length(Tq_conv)*Ts;
PxXi_mod = (abs(fftshift(fft(Xi_mod,N))).^2)./length(Ti_conv)*Ts;

figure;
subplot(2,1,1)
plot(F,PxXq_mod);
title('Periodogram of PxXq_mod');
xlabel('F(Hz)');
ylabel('PxXq_mod');

subplot(2,1,2);
plot(F,PxXi_mod);
title('Periodogram of PxXi_mod');
xlabel('F(Hz)');
ylabel('PxXi_mod');

figure;
subplot(2,1,1)
semilogy(F,PxXq_mod);
title('Periodogram of PxXq_mod with semilogy');
xlabel('F(Hz)');
ylabel('PxXq_mod');

subplot(2,1,2);
semilogy(F,PxXi_mod);
title('Periodogram of PxXi_mod with semilogy');
xlabel('F(Hz)');

%a.6
X_mod = Xi_mod + Xq_mod;

figure;
plot(Tq_conv,X_mod);
title('X_mod');
xlabel('time');
ylabel('X_mod');

PxX_mod = (abs(fftshift(fft(X_mod,N))).^2)./length(Tq_conv)*Ts;

%periodogram X_mod
figure;
subplot(2,1,1);
plot(F,PxX_mod);
title('Periodogram of PxX_mod');
xlabel('F(Hz)');
ylabel('PxX_mod');

subplot(2,1,2);
semilogy(F,PxX_mod);
title('Periodogram of PxX_mod using semilogy');
xlabel('F(Hz)');
ylabel('PxX_mod');


%a.8
SNR = 20;
varW = (10*A^2)/(Ts*(10^(SNR/10)));
W_G_signal = sqrt(varW)*randn(1,length(X_mod));
W = W_G_signal + X_mod;

%a.9
Wi = W.*cos(2*pi*F0*Ti_conv);
Wq = W.*(-sin(2*pi*F0*Tq_conv));

figure;
subplot(2,1,1);
plot(Ti_conv ,Wi);
title('Wi with white noise');
ylabel('Wi');
xlabel('time');

subplot(2,1,2);
plot(Tq_conv ,Wq);
title('Wq with white noise');
ylabel('Wq');
xlabel('time');

PxWi = (abs(fftshift(fft(Wi,N))).^2)./length(Ti_conv)*Ts;
PxWq = (abs(fftshift(fft(Wq,N))).^2)./length(Tq_conv)*Ts;

figure;
subplot(2,1,1);
plot(F,PxWi);
title('Periodogram of PxWi');
xlabel('F(Hz)');
ylabel('PxWi');

subplot(2,1,2);
plot(F,PxWq);
title('Periodogram of PxWq');
xlabel('F(Hz)');
ylabel('PxWq');


figure;
subplot(2,1,1);
semilogy(F,PxWi);
title('Periodogram of PxWi using semilogy');
xlabel('F(Hz)');
ylabel('PxWi');

subplot(2,1,2);
semilogy(F,PxWq);
title('Periodogram of PxWq using semilogy');
xlabel('F(Hz)');
ylabel('PxWq');


%a.10
Wi_delta = conv(Wi,phi)*Ts;
Wq_delta = conv(Wq,phi)*Ts;

Ti_conv_new = linspace(Ti_conv(1)+t(1),Ti_conv(end)+t(end),length(Wi_delta));
Tq_conv_new = linspace(Tq_conv(1)+t(1),Tq_conv(end)+t(end),length(Wq_delta));

figure;
subplot(2,1,1);
plot(Ti_conv_new,Wi_delta);
title('Convolution of phi with Wi');
xlabel('time');
ylabel('Wi_delta');

subplot(2,1,2);
plot(Tq_conv_new,Wq_delta);
title('Convolution of phi with Wq');
xlabel('time');
ylabel('Wq_delta');

PxWi_new = (abs(fftshift(fft(Wi_delta,N))).^2)./length(Ti_conv_new)*Ts;
PxWq_new = (abs(fftshift(fft(Wq_delta,N))).^2)./length(Tq_conv_new)*Ts;

figure;
subplot(2,1,1);
plot(F,PxWi_new);
title('Periodogram of conv between phi and Wi');
xlabel('F(Hz)');
ylabel('PxWi_new');

subplot(2,1,2);
plot(F,PxWq_new);
title('Periodogram of conv between phi and Wq');
xlabel('F(Hz)');
ylabel('PxWq_new');


figure;
subplot(2,1,1);
semilogy(F,PxWi_new);
title('Periodogram of conv between phi and Wi using semilogy');
xlabel('F(Hz)');
ylabel('PxWi_new');

subplot(2,1,2);
semilogy(F,PxWq_new);
title('Periodogram of conv between phi and Wq using semilogy');
xlabel('F(Hz)');
ylabel('PxWq_new');

%a.11
Zi_downsampled = downsample(Wi_delta,over);
Zq_downsampled = downsample(Wq_delta,over);

figure();
scatter(Zi_downsampled,Zq_downsampled);
title('Symbols of Z');
ylabel('Z(t)')

%a.12
Yi_est = detect_4_PAM(Zi_downsampled,A);
Yq_est = detect_4_PAM(Zq_downsampled,A);

figure;
scatter(Yi_est,Yq_est);
title('Symbols of Y');

%a.13
counterq=0;
counteri=0
for i=1:length(Xi)
%     if(Yi_est ~= Yq_est)
    if(Yi_est(i) ~= round(Xi(i)))     
        counteri = counteri +1;
    end
    if(Yq_est(i) ~= round(Xq(i)))     
        counterq = counterq +1;
    end
end
num_of_errors = counterq+counteri

%a.14
i_est_bit = PAM_4_to_bits(Yi_est,A);
q_est_bit = PAM_4_to_bits(Yq_est,A);

N=500;

est_bit = zeros(1,4*(N+16));
est_bit(1:2*(N+16)) = i_est_bit;
est_bit(2*N+32+1:end) = q_est_bit;

%a.15
error = 0;

for i=1:length(b)
    if(b(i) ~= est_bit(i))
        error = error +1;
    end
end
num_of_errors = error


%b.1
SNRdb = 0:2:16;
K = 200;

for n=1:length(SNRdb)
    num_of_error = 0;
    bit_error = 0;
    for j=1:K
        varW = (10*A^2)/(Ts*(10^(SNRdb(n)/10)));
        W_G_signal = sqrt(varW)*randn(1,length(X_mod));
        W = W_G_signal + X_mod;
        
        Wi = W.*cos(2*pi*F0*Ti_conv);
        Wq = W.*(-sin(2*pi*F0*Tq_conv));
        
        Wi_delta = conv(Wi,phi)*Ts;
        Wq_delta = conv(Wq,phi)*Ts;
        Zi_downsampled = downsample(Wi_delta,over);
        Zq_downsampled = downsample(Wq_delta,over);
        
        Yi_est = detect_4_PAM(Zi_downsampled,A);
        Yq_est = detect_4_PAM(Zq_downsampled,A);
        
        counterq = 0;
        counteri = 0;
        for i=1:length(Xi)
        %     if(Yi_est ~= Yq_est)
            if(Yi_est(i) ~= round(Xi(i)))     
                counteri = counteri +1;
            end
            if(Yq_est(i) ~= round(Xq(i)))     
                counterq = counterq +1;
            end
        end
        num_of_errors = num_of_errors + counterq + counteri;

        i_est_bit = PAM_4_to_bits(Yi_est,A);
        q_est_bit = PAM_4_to_bits(Yq_est,A);
        
        est_bit = zeros(1,4*(N+16));
        est_bit(1:2*(N+16)) = i_est_bit;
        est_bit(2*N+32+1:end) = q_est_bit;
        
        for i=1:length(b)
            if(b(i) ~= est_bit(i))
                bit_error = bit_error + 1;
            end
        end
    end
    
    err_exp(1,n) = num_of_error/(N*K);
    ber_exp(1,n) = bit_error/(N*K*4);
end

Theo_symbol = 3/2.*erfc(sqrt(0.1*(10.^(SNRdb/10))))-(1/4)*3/2.*erfc(sqrt(0.1*(10.^(SNRdb/10))));
Theo_bit = (1/4)*3/2.*Theo_symbol;

%b.2
figure;
semilogy(SNRdb,err_exp);
title('SER using Monte Carlo Method');
xlabel('SNRdb');
ylabel('SER for 16-QAM');

hold on;
semilogy(SNRdb,Theo_symbol);

legend('Experimental','Theoritical');

%b.3
figure;
semilogy(SNRdb,ber_exp);
title('BER using Monte Carlo Method');
xlabel('SNRdb');
ylabel('BER for 16-QAM');

hold on;
semilogy(SNRdb,Theo_bit);

legend('Experimental','Theoritical');
        
