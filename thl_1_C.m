%c
clc;
close all;
clear all;

%c1
a = 0.5;
A = 5;
T = 0.1;
over = 10;
Ts = T/over;
N = 100;
 
%creation of the N=50 random bits
b = (sign(randn(N, 1)) + 1)/2; 

%c2
%a
%transforming the N bits in to 2-pam
X = bits_to_2PAM(b);

%b
%simulation of the signal X{delta}(t)
X_delta = 1/Ts*upsample(X, over);
 
%time
t = 0:Ts:((N+N*(over-1))-1)*Ts; %exoume prosthiki over-1 midenikwn meta3i simvolwn
 
figure;
plot(t,X_delta);
xlabel('t');
title('X{delta}(t) = sum(Xk*delta(t-kT))');
 
%creation of fi(t) pulse
[phi_1,t1] = srrc_pulse(T, Ts, A, a);

%time of convolution
tx = t(1)+t1(1) :Ts: t(end)+t1(end);

%conv function for the convolution between fi(t) and X_delta
X_conv = conv(X_delta,phi_1)*Ts;
 
figure;
plot(tx,X_conv);
xlabel('t');
title('X(t) = X{delta}(t) * fi(t)')
 
%creation of fi(-t) pulse
phi_2 = phi_1(end:-1:1);
t2 = -t1(end:-1:1);

%convolution of fi(-t) and X_conv
z = conv(X_conv,phi_2)*Ts;

%time of the convolution
tz= tx(1)+t2(1) :Ts: tx(end)+t2(end);
 
figure;
plot(tz,z);
xlabel('t');
title('Z(t)=X(t)*fi(-t)')
 
figure;
plot(tz,z);
hold on;
%plot in the same figure Z and X
stem([0:N-1]*T,X,'r');
xlabel('t');
legend('Z(t)=X(t)*fi(-t)','X{k}')
 