clc;
close all;
clear all;

T = 0.001;
over = 10;
Ts = T/over;
Fs = 1/Ts;
A = 4;
a = 0.5;
Nf = 2048;
N = 100;
f0 = 3000; %frequency
 
%srrc pulse)
[phi, t] = srrc_pulse(T, Ts, A, a); 
 
%fourier transforamtion
y=abs(fftshift(fft(phi,Nf)*Ts));
 
%frequency for the bandwidth
F = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;
 
%creation of the bits
b = (sign(randn(N, 1)) + 1)/2;
 
%2-PAM modulation
X=bits_to_2PAM(b);
delta_X = 1/Ts*upsample(X, over);
 
%time of the conv 
t1 = 0:Ts:((N+N*(over-1))-1)*Ts;

%convolution
x = conv(delta_X,phi)*Ts;
tx = t(1)+t1(1) :Ts: t(end) + t1(end); 
%Theta uniform in [0,2pi)
theta = 2*pi*rand(1);
 
%random signal Y(t)
Yt = x.*(cos(2*pi*f0.*transpose(tx) + theta));
 
Py = (abs(fftshift(fft(Yt,Nf)*Ts)).^2)/(length(Yt)*Ts); %pasmatikh piknothta isxios
 
figure;
subplot(2,1,1);
plot(F,Py);
xlabel('Frequency(Hz)');
ylabel('Py(F)');
title('Periodic diagram');
 
subplot(2,1,2);
semilogy(F,Py);
xlabel('Frequency (Hz)');
xlabel('Py(F)');
title('Periodic diagram');

K=1000;
Py_new = zeros(K,Nf);
for i=1:K 
    b = (sign(randn(N, 1)) + 1)/2;    
    X = bits_to_2PAM(b);
    delta_X = 1/Ts*upsample(X, over); 
    tx = t(1)+t1(1) :Ts: t(end) + t1(end); 
    x = conv(delta_X,phi)*Ts;
    theta = 2*pi*rand(1);
    yt = x.*(cos(2*pi*f0.*transpose(tx)+theta));
    Py_new(i,:) = (abs(fftshift(fft(yt,Nf)*Ts)).^2)/(length(yt)*Ts);
end;
 
%fasmatikh piknothta isxios
matlab=sum(Py_new)./K;
 
figure;
semilogy(F,matlab)
xlabel('Frequency(Hz)');
ylabel('Py');
title('Fasmatikh piknothta isxios');
legend('matlab');
