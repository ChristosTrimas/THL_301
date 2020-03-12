clc;
close all;
clear all;

%a.1

T = 10^(-3);
over = 10;
Ts = T / over;
A = 4;
a = 0.5;
Fs = 1/Ts;
Nf = 2048; %from the exercise

%SRRC from lab 1
[phi, t] = srrc_pulse(T, Ts, A, a); 
 
%Fourier Transformation 
y=abs(fftshift(fft(phi,Nf)*Ts)); %fftshift centeralizes the transformation in zero
 
f = -Fs/2:Fs/Nf:Fs/2-Fs/Nf; %frequency like lab 1
 
%Ipologismos famsitkhs piknothtas
figure;
semilogy(f,y.^2,'b');
xlabel('Frequency (Hz)');
ylabel('|FI(F)|');
title('Semilogy famstikh piknothta energeias');
legend('|FI(F)|^2');

%a.2

N = 100;

%creation of the bits, independent and with the same probability
b = (sign(randn(N, 1)) + 1)/2;
 
%Use of the function that was created in the 1st exercise to make the
%mapping
X = bits_to_2PAM(b);
 
%like the C section of lab1
X_delta = 1/Ts*upsample(X, over);
 
%time 
t1 = 0:Ts:((N+N*(over-1))-1)*Ts;
 
%Convolutionof X_delta fi(t)
x = conv(X_delta,phi)*Ts;

%time of convolution 
tc = t(1)+t1(1) :Ts: t(end) + t1(end);

figure;
plot(tc,x);
xlabel('t');
title('X(t)');

%a.3
Px=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tc)*Ts);
 
figure;
subplot(2,1,1);
plot(f,Px)
xlabel('Frequency(Hz)');
ylabel('X(t)');
title('Periodic diagram of X(t)');
 
subplot(2,1,2);
semilogy(f,Px) %logarithmimeno
xlabel('Frequency(Hz)');
ylabel('X(t)');
title('Periodic diagram of X(t)');

K = 1000;
Pxx = zeros(K,Nf);

for i = 1:K 
    b = (sign(randn(N, 1)) + 1)/2; %creation of the bits like before
    %2 pam modulation
    X = bits_to_2PAM(b);
    X_delta = 1/Ts*upsample(X, over); %like before
    %time 
    tx = 0:Ts:((N+N*(over-1))-1)*Ts;
    %convolution
    x = conv(X_delta,phi)*Ts;
 
    Pxx(i,:)=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
end;
 
matlab_Px = sum(Pxx) ./ K;
theoritical_Px = (var(X)^2/T).*(y.^2); %theoritical
 
figure;
semilogy(f,matlab_Px,'b')
hold on;
semilogy(f,theoritical_Px,'r')
xlabel('Frequency (Hz)');
title('Fasmatikh piknothta isxuos,K=1000');
legend('matlab','theoritical');
 
%a.4

b = (sign(randn(N, 1)) + 1)/2; %like before
 
X = bits_to_4PAM(b); %use of the created function bits_to_4PAM
 
%like lab1
X_delta = 1/Ts*upsample(X, over);
 
tx = 0:Ts:((N+N*(over-1))-1)*Ts;
 
%convolution of X_delta
x = conv(X_delta,phi)*Ts;
 
Px =(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);

figure;
subplot(2,1,1);
plot(f,Px)
xlabel('Frequency(Hz)');
ylabel('X(t)');
title('Periodic diagram of X(t)');
 
subplot(2,1,2);
semilogy(f,Px) %logarithmimeno
xlabel('Frequency(Hz)');
ylabel('X(t)');
title('Periodic diagram of X(t)');

Pxxx=zeros(K,Nf);

for i=1:K 
    
    b = (sign(randn(N, 1)) + 1)/2;
    X = bits_to_4PAM(b);
    X_delta = 1/Ts*upsample(X, over);
    tx = 0:Ts:((N+N*(over-1))-1)*Ts;
    x = conv(X_delta,phi)*Ts;
    Pxxx(i,:) = (abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
end;

theoritical = (5/T).*(y.^2);
matlab = sum(Pxxx)./K;
 
figure;
semilogy(f,matlab,'b')
hold on;
semilogy(f,theoritical,'r')
xlabel('Frequency (Hz)');
title('Fasmatikh Piknothta Isxios semilogy');
legend('matlab','theoritical');

%a.5
T_new = 2*T;
over_new = 2*over;
Ts_new = T_new / over_new; %same Ts with a.2
Fs_new = 1 / Ts_new;

%new srrc pulse
[phi,t] = srrc_pulse(T, Ts, A, a);

z = abs(fftshift(fft(phi,Nf)*Ts_new));

%frequency vector for the whole bandwidth
F = -Fs_new/2 :Fs_new/Nf: Fs_new/2 - Fs/Nf;

b = (sign(randn(N,1)) + 1)/2;
X = bits_to_2PAM(b);
delta_X = 1/Ts*upsample(X,over_new);
tx = 0:Ts_new:((N+N*(over_new-1))-1)*Ts_new;
x = conv(delta_X,phi)*Ts;
Px_new = (abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);

figure;
subplot(2,1,1);
plot(F,Px_new);
xlabel('Frequency(Hz)');
ylabel('Px');
title('Periodic Diagram of X');

subplot(2,1,2);
semilogy(F,Px_new);
xlabel('Frequency(Hz)');
ylabel('Px');
title('Periodic diagram of X(t)');

Px = zeros(K,Nf);
% K=100;
% N=50;
for i = 1:K
    b = (sign(randn (N,1)) + 1)/2;
    X = bits_to_2PAM(b);
    delta_X = 1/Ts_new*upsample(X,over_new);
    tx = 0:Ts_new:((N+N*(over_new-1))-1)*Ts_new; 
    x = conv(delta_X,phi)*Ts_new;
    Px(i,:)=(abs(fftshift(fft(x,Nf)*Ts_new)).^2)/(length(tx)*Ts_new);
end;

matlab = sum(Px) ./ K;
theoritical = (var(X)^2 / T) .* (z.^2);

figure;
semilogy(F,matlab,'b');
hold on;
semilogy(F,theoritical,'r');
xlabel('Frequency(Hz)');
ylabel('Px');
title('Fasmatikh piknothta isxios');
legend('matlab','theoritical');