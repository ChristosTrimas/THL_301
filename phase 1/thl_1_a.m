clc;
close all;
clear all;
%a1

T = 10^(-2); %period
over = 10;
Ts = T/over; %periodos deigmatolipsias
A = 4;
a1 = 0;   %green pulse
a2 = 0.5; %blue pulse
a3 = 1;   %red pulse
 
 
%creation of srrc pulses
[phi1,t1] = srrc_pulse(T, Ts, A, a1); 
[phi2,t2] = srrc_pulse(T, Ts, A, a2); 
[phi3,t3] = srrc_pulse(T, Ts, A, a3);
 
%plots
figure;    
plot(t1,phi1,'g')
hold on
plot(t2,phi2,'b')
hold on
plot(t3,phi3,'r')
xlabel('t')
ylabel('fi(t)')
title('3 different srrc pulses')
legend('a = 0','a = 0.5', 'a =1')
grid on 

%a2
Fs = 1/Ts;    
Nf = 2048;                      %samples
F = [-Fs/2:Fs/Nf:Fs/2-Fs/Nf]; %frequency
 
y1 = fftshift(fft(phi1,Nf)*Ts); %fourier transformation
y2 = fftshift(fft(phi2,Nf)*Ts);
y3 = fftshift(fft(phi3,Nf)*Ts);
 
figure;
plot(F,abs(y1).^2,'g')
hold on
plot(F,abs(y2).^2,'b')
hold on
plot(F,abs(y3).^2,'r')
xlabel('F');
ylabel('|FI(F)|^2')
title('Fasmatikh Piknothta Isxios')
legend('a = 0','a = 0.5', 'a =1')
 
figure;
semilogy(F,abs(y1).^2,'g');
hold on;
semilogy(F,abs(y2).^2,'b');
hold on;
semilogy(F,abs(y3).^2,'r');
xlabel('F');
ylabel('|FI(F)|^2')
title('Logarithmimenh Fasmatikh Piknothta Isxuos')
legend('a = 0','a = 0.5', 'a =1')
 
%a3
 
%bandwith based in theory
BW1 = (1+a1)./(2*T)
BW2 = (1+a2)./(2*T)
BW3 = (1+a3)./(2*T)
 
%c1,c2
c1 = T./(10^3);
c2 = T./(10^5);
 
%emfanisi fasmatikis piknotias energeias(semilogy) mazi me tin c1
figure;
semilogy(F,abs(y1).^2,'g');
hold on;
semilogy(F,abs(y2).^2,'b');
hold on;
semilogy(F,abs(y3).^2,'r');
hold on;
plot(F, c1, 'k.')
hold on;
xlabel('F(Hz)');
ylabel('FI(F)^2')
legend('a = 0','a = 0.5', 'a =1', 'c=T/10^3')
 
%emfanisi fasmatikis piknotias energeias(semilogy) mazi me tin c2
figure;
semilogy(F,abs(y1).^2,'g');
hold on;
semilogy(F,abs(y2).^2,'b');
hold on;
semilogy(F,abs(y3).^2,'r');
hold on;
plot(F, c2, 'c.')
hold on;
xlabel('F(Hz)');
ylabel('FI(F)^2')
legend('a = 0','a = 0.5', 'a =1','c=T/10^5') 
