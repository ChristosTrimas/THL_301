%b
clc;
clear all;
close all;

T = 10^(-2);
over = 10;
Ts = T/over;
A = 5;
a1= 0;
a2 = 0.5;
a3 = 1;
 
 
%dimiourgia palmwn srrc
[phi_1,t1] = srrc_pulse(T, Ts, A, a1);
[phi_2,t2] = srrc_pulse(T, Ts, A, a2);
[phi_3,t3] = srrc_pulse(T, Ts, A, a3);
 
%a=0
for k = setdiff(0:4,3)
    %zero padding to create the signal moved by kT
    phi_1_kt = [zeros(1,(1/Ts)*k*T) phi_1(1:end-(1/Ts)*k*T)]; 
    z1 = phi_1.*phi_1_kt; %ginomeno 
    
    figure;
    subplot(2,1,1);
    plot(t1,phi_1,'b');
    hold on;
    plot(t1,phi_1_kt,'r');
    grid on;
    xlabel('t');
    
    if(k==0)
        title('fi(t), fi(t-kT) for a=0 and k=0')
    elseif(k==1)
        title('fi(t), fi(t-kT) for a=0 and k=1')
    elseif(k==2)
        title('fi(t), fi(t-kT) for a=0 and k=2')
    elseif(k==4)
        title('fi(t), fi(t-kT) for a=0 and k=4')
    end;
    legend('fi(t)', 'fi(t-kT)')   
    
    subplot(2,1,2);
    plot(t1,z1); %ginomeno
    xlabel('t');
    
    if(k==0)
        title('fi(t), fi(t-kT) for a=0 and k=0')
    elseif(k==1)
        title('fi(t), fi(t-kT) for a=0 and k=1')
    elseif(k==2)
        title('fi(t), fi(t-kT) for a=0 and k=2')
    elseif(k==4)
        title('fi(t), fi(t-kT) for a=0 and k=4')
    end;
    legend('fi(t) X fi(t-kT)');
    
    integral_1 = sum(z1)*Ts %ypologismos oloklirwmatos
    grid on
    
end;
 
%a=0,5
for k = setdiff(0:4,3)
    phi_2_kt = [zeros(1,(1/Ts)*k*T) phi_2(1:end-(1/Ts)*k*T)];
    z2 = phi_2.*phi_2_kt;
    
    figure;
    subplot(2,1,1)
    plot(t2,phi_2,'b')
    hold on;
    plot(t2,phi_2_kt,'r')
    grid on
    xlabel('t')
    
    if(k==0)
        title('fi(t), fi(t-kT) for a=0.5 and k=0')
    elseif(k==1)
        title('fi(t), fi(t-kT) for a=0.5 and k=1')
    elseif(k==2)
        title('fi(t), fi(t-kT) for a=0.5 and k=2')
    elseif(k==4)
        title('fi(t), fi(t-kT) for a=0.5 and k=4')
    end;
    legend('fi(t)', 'fi(t-kT)')
    
    subplot(2,1,2)
    plot(t2,z2)
    xlabel('t')
    
    if(k==0)
        title('fi(t), fi(t-kT) for a=0.5 and k=0')
    elseif(k==1)
        title('fi(t), fi(t-kT) for a=0.5 and k=1')
    elseif(k==2)
        title('fi(t), fi(t-kT) for a=0.5 and k=2')
    elseif(k==4)
        title('fi(t), fi(t-kT) for a=0.5 and k=4')
    end;
    legend('fi(t) X fi(t-kT)')
    
    integral_2 = sum(z2)*Ts
    grid on
    
end;
 
%a=1
for k = setdiff(0:4,3)
    phi_3_kt = [zeros(1,(1/Ts)*k*T) phi_3(1:end-(1/Ts)*k*T)];
    z3 = phi_3.*phi_3_kt;
    
    figure;
    subplot(2,1,1);
    plot(t3,phi_3,'b');
    hold on;
    plot(t3,phi_3_kt,'r');
    grid on;
    xlabel('t');
    
    if(k==0)
        title('fi(t), fi(t-kT) for a=1 and k=0')
    elseif(k==1)
        title('fi(t), fi(t-kT) for a=1 and k=1')
    elseif(k==2)
        title('fi(t), fi(t-kT) for a=1 and k=2')
    elseif(k==4)
        title('fi(t), fi(t-kT) for a=1 and k=4')
    end;
    legend('fi(t)', 'fi(t-kT)')
    
    subplot(2,1,2);
    plot(t3,z3);
    xlabel('t');
    
    if(k==0)
        title('fi(t), fi(t-kT) for a=1 and k=0')
    elseif(k==1)
        title('fi(t), fi(t-kT) for a=1 and k=1')
    elseif(k==2)
        title('fi(t), fi(t-kT) for a=1 and k=2')
    elseif(k==4)
        title('fi(t), fi(t-kT) for a=1 and k=4')
    end;
    legend('fi(t) X fi(t-kT)')
    
    integral_3 = sum(z3)*Ts
    grid on;
    
end;
 
