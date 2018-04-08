%% Assignment 4 - Kyle Poulin 100939284
%% Question 1
% 1:
%
% Time & Frequency: $V_i = V_{in}$
%
% 2:
%
% Time: $G_1(V_2-V_1)+C(\frac{d(V_2-V_1)}{dt} )+G_2V_2-I_L=0$
%
% Frequency:    $G_1(V_2-V_1)+C(j \omega (V_2-V_1) )+G_2V_2-I_L=0$
%
% 3:
%
% Time: $$V_2-V_3-L \frac {dI_l}{dt}=0$$
%
% Frequency:    $$V_2-V_3-L (j \omega)I_l=0$$
%
% 4:
%
% Time & Fequency:  $$-I_L+G_3V_3=0$$
%
% 5:
%
% Time & Fequency:  $$V_4-\alpha I_3=0$$
%
% 6:
%
% Time & Fequency:  $$G_3V_3-I_3=0$$
%
% 7:
%
% Time & Fequency:  $$G_4(V_0-V_4)+G_0 V_0 =0$$
%


C = 0.25;
Cconst = 0.25;
Cstd = 0.05;
L = 0.2;
G1 = 1;
G2 = 0.5;
G3 = 0.1;
G4 = 10;
G0 = 1/1000;
a = 100;
Vin=1;


Cm = [0 0 0 0 0 0 0; -C C 0 0 0 0 0; 0 0 -L 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0];
Gm = [1 0 0 0 0 0 0; -G2 G1+G2 -1 0 0 0 0; 0 1 0 -1 0 0 0; 0 0 -1 G3 0 0 0; 0 0 0 0 -a 1 0; 0 0 0 G3 -1 0 0; 0 0 0 0 0 -G4 G4+G0];
%V = [V1; V2; IL; V3; I3; V4; V0];
F = [Vin; 0; 0; 0; 0; 0; 0];

V0 = zeros(1,100);
V3 = zeros(1,100);
counter = 1;

for i= -10:0.2:10
   
   F = [i; 0; 0; 0; 0; 0; 0];
   V = Gm\F;
   V0(1,counter) = V(7);
   V3(1,counter) = V(1);
   counter = counter+1;
end

figure(1);
plot(V3);
title('V3 vs Vin for DC 0-10 V input');
xlabel('Input (V)');
ylabel('V3 (V)');
figure(2);
plot(V0);
title('Ouput vs. Vin for DC 0-10 V input')
xlabel('Input (V)')
ylabel('Output (V)')

F = [Vin; 0; 0; 0; 0; 0; 0];


xaxis = linspace(0,1000,1001);
gain = [];
for i = 0:1000
    V = (Gm+1j*i*Cm)\F;
    gain = [gain, V(7)];
end

%plot(gain);
figure(3);
semilogx(xaxis,gain);
title('Output vs. frequency')
xlabel('freq(rad/s)')
ylabel('Ouput (V)')
figure(4);
semilogx(xaxis,20*log10(gain));
title('Output vs. frequency')
xlabel('freq(rad/s)')
ylabel('Ouput (V)')

F = [Vin; 0; 0; 0; 0; 0; 0];
Gm = [1 0 0 0 0 0 0; -G2 G1+G2 -1 0 0 0 0; 0 1 0 -1 0 0 0; 0 0 -1 G3 0 0 0; 0 0 0 0 -a 1 0; 0 0 0 G3 -1 0 0; 0 0 0 0 0 -G4 G4+G0];
gain = [];
for i=1:10000
    C = normrnd(Cconst,Cstd);
    Cm = [0 0 0 0 0 0 0; -C C 0 0 0 0 0; 0 0 -L 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0];
    V = (Gm+pi*Cm)\F;
    gain = [gain, V(7)];
end

figure(5);
histogram(real(20*log10(gain)));
title('Gain Hist')
xlabel('Gain (dB)')
ylabel('Counts')

%% Question 2


% Step Response:
vin = 1;
F1 = 0;
V1=0;
t = linspace(0,1,1000);
dt = 0.001;


for i=1:31
    F1(i,1:7) = [0;0;0;0;0;0;0];
end
for i=32:1000
    F1(i,1:7) = [vin;0;0;0;0;0;0];
end

V1(1:7,1) = (Cm/dt+Gm)^-1 *(F1(1,:)');
for i=2:1000
    V1(:,i) = (Cm/dt+Gm)^-1 *(Cm*V1(:,i-1)/dt+F1(i,:)');
end

figure(6);
plot(t,V1(7,:));
hold on
plot(t,F1(:,1));
title('Step Response')
xlabel('t (s)')
ylabel('Voltage (V)')
legend('Output','Voltage Step')

figure(7)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(V1(7,:)))))
hold on 
semilogy(linspace(-500,500,1000),fftshift(abs(fft(F1(:,1)))))

title('Step Response Fourier Transform')
xlabel('Freq (rad/s)')
ylabel('Voltage (V)')
legend('Output','input signal')

% Sine input
V2 = 0;
freq = 1/0.03;
F2=0;
for i=1:1000
    F2(i,1:7) = [sin(2*pi*freq*t(i)),0,0,0,0,0,0];
end


V2(1:7,1) = (Cm/dt+Gm)^-1 *(F2(1,:)');
for i=2:1000
    V2(:,i)=(Cm/dt+Gm)^-1 *(Cm*V2(:,i-1)/dt+F2(i,:)');   
end

% Plot the response
figure(8)
plot(t,V2(7,:))
hold on 
plot(t,F2(:,1))

title('Sine Response')
xlabel('t (s)')
ylabel('Voltage (V)')
legend('Output','Input')

figure(9)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(V2(7,:)))))
hold on 
semilogy(linspace(-500,500,1000),fftshift(abs(fft(F2(:,1)))))

title('Sine Response Fourier Transform')
xlabel('Freq (rad/s)')
ylabel('Voltage (V)')
legend('Output','input signal')

% Gauss response
V3 = 0;
F3 = 0;
for i=1:1000
    F3(i,1:7)=[exp(-1/2*((t(i)-0.06)*freq)^2),0,0,0,0,0,0];
end

V3(1:7,1) = (Cm/dt+Gm)^-1 *(F3(1,:)');
for i=2:1000
    V3(:,i)=(Cm/dt+Gm)^-1 *(Cm*V3(:,i-1)/dt+F3(i,:)');   
end


figure(10)
plot(t,V3(7,:))
hold on 
plot(t,F3(:,1))

title('Guassian Response')
xlabel('t (s)')
ylabel('Voltage (V)')
legend('Output','Input')

figure(11)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(V3(7,:)))))
hold on 
semilogy(linspace(-500,500,1000),fftshift(abs(fft(F3(:,1)))))

title('Guassian Response Fourier Transform')
xlabel('Freq (rad/s)')
ylabel('Voltage (V)')
legend('Output','input signal')

% Larger Timestep
V4 = 0;
t = linspace(0,1,100);

V4(1:7,1)=(Cm/dt+Gm)^-1 *(F3(1,:)');
for i=2:100
    V4(:,i)=(Cm/dt+Gm)^-1 *(Cm*V3(:,i-1)/dt+F3(i,:)');   
end

figure(12)
plot(t,V4(7,:))
hold on 
t = linspace(0,1,1000);
plot(t,F3(:,1))
legend('Vout','input signal')
title('DC Guassian responce')
xlabel('Time (s)')
ylabel('Voltage (V)')

% The longer time step can be seen in figure 12. When acting on the
% gaussian input, the peak widens and gets pushed back in time.

%% Question 3


% Redefine C and G matrixes to be 8x8 to account for the new junction
In = 0.001*randn();
Cn = 0.00001;

Cm = [0 0 0 0 0 0 0 0; -C C 0 0 0 0 0 0; 0 0 -L 0 0 0 0 0; 0 0 0 Cn 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0];
Gm = [1 0 0 0 0 0 0 0; -G2 G1+G2 -1 0 0 0 0 0; 0 1 0 -1 0 0 0 0; 0 0 -1 G3 0 0 0 -1; 0 0 0 0 -a 1 0 0; 0 0 0 G3 -1 0 0 0; 0 0 0 0 0 -G4 G4+G0 0; 0 0 0 0 0 0 0 1];
F = [Vin; 0; 0; 0; 0; 0; 0; In];

V1 = 0;
F1 = 0;

for i=1:1000
    In = dt*randn();
    F1(i,1:8)=[exp(-1/2*((t(i)-0.06)*freq)^2),0,0,0,0,0,0,In];
end

V1(1:8,1)=(Cm/dt+Gm)^-1*(F1(1,:)');
for i=2:1000
    V1(:,i)=(Cm/dt+Gm)^-1*(Cm*V1(:,i-1)/dt+F1(i,:)');   
end

figure(13)
plot(t,V1(7,:))
hold on 
plot(t,F1(:,1))

title('Gaussian with Noise')
xlabel('t (s)')
ylabel('Voltage (V)')
legend('Output','Input')

figure(14)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(V1(7,:)))))
hold on 
semilogy(linspace(-500,500,1000),fftshift(abs(fft(F3(:,1)))))

title('Guassian with Noise Fourier Transform')
xlabel('Freq (rad/s)')
ylabel('Voltage (V)')
legend('Output','Input')


Cn = 0.00002;
Cm = [0 0 0 0 0 0 0 0; -C C 0 0 0 0 0 0; 0 0 -L 0 0 0 0 0; 0 0 0 Cn 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0];

V2(1:8,1)=(Cm/dt+Gm)^-1*(F1(1,:)');
for i=2:1000
    V2(:,i)=(Cm/dt+Gm)^-1*(Cm*V2(:,i-1)/dt+F1(i,:)');   
end

figure(15)
plot(t,V2(7,:))
hold on 
plot(t,F1(:,1))

title('Gaussian with Noise (Cn = 20uF)')
xlabel('t (s)')
ylabel('Voltage (V)')
legend('Output','Input')


Cn = 0.0002;
Cm = [0 0 0 0 0 0 0 0; -C C 0 0 0 0 0 0; 0 0 -L 0 0 0 0 0; 0 0 0 Cn 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0];

V3(1:8,1)=(Cm/dt+Gm)^-1*(F1(1,:)');
for i=2:1000
    V3(:,i)=(Cm/dt+Gm)^-1*(Cm*V3(:,i-1)/dt+F1(i,:)');   
end

figure(16)
plot(t,V3(7,:))
hold on 
plot(t,F1(:,1))

title('Gaussian with Noise (Cn = 200uF)')
xlabel('t (s)')
ylabel('Voltage (V)')
legend('Output','Input')


Cn = 0.002;
Cm = [0 0 0 0 0 0 0 0; -C C 0 0 0 0 0 0; 0 0 -L 0 0 0 0 0; 0 0 0 Cn 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0];

V4(1:8,1)=(Cm/dt+Gm)^-1*(F1(1,:)');
for i=2:1000
    V4(:,i)=(Cm/dt+Gm)^-1*(Cm*V4(:,i-1)/dt+F1(i,:)');   
end

figure(17)
plot(t,V4(7,:))
hold on 
plot(t,F1(:,1))

title('Gaussian with Noise (Cn = 2mF)')
xlabel('t (s)')
ylabel('Voltage (V)')
legend('Output','Input')

% Return the Cn value to default
Cn = 0.00001;
Cm = [0 0 0 0 0 0 0 0; -C C 0 0 0 0 0 0; 0 0 -L 0 0 0 0 0; 0 0 0 Cn 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0];

% It can be seen that increasing the capacitor value decreases the noise.
% Once the value has been increased to 2mF, the noise is negligable.

V5 = 0;
F2 = 0;

t = linspace(0,1,800);
dt = 0.00125;

for i=1:800
    In = 0.001*randn();
    F2(i,1:8)=[exp(-1/2*((t(i)-0.06)*freq)^2),0,0,0,0,0,0,In];
end

V5(1:8,1)=(Cm/dt+Gm)^-1*(F2(1,:)');
for i=2:800
    V5(:,i)=(Cm/dt+Gm)^-1*(Cm*V5(:,i-1)/dt+F2(i,:)');   
end

figure(18)
plot(t,V5(7,:))
hold on 
plot(t,F2(:,1))

title('Gaussian with Noise (dt = 1.25ms)')
xlabel('t (s)')
ylabel('Voltage (V)')
legend('Output','Input')



V6 = 0;
F3 = 0;

t = linspace(0,1,500);
dt = 0.1;

for i=1:500
    In = 0.001*randn();
    F3(i,1:8)=[exp(-1/2*((t(i)-0.06)*freq)^2),0,0,0,0,0,0,In];
end

V6(1:8,1)=(Cm/dt+Gm)^-1*(F3(1,:)');
for i=2:500
    V6(:,i)=(Cm/dt+Gm)^-1*(Cm*V6(:,i-1)/dt+F3(i,:)');   
end

figure(19)
plot(t,V6(7,:))
hold on 
plot(t,F3(:,1))

title('Gaussian with Noise (dt = 0.1s)')
xlabel('t (s)')
ylabel('Voltage (V)')
legend('Output','Input')

% With an increase in the timestep, we see that the output shape does not
% dip as far into the negative voltage range. The output signal is also
% slightly less noisy (figure 19) as compared to that of a smaller time
% step (figure 18).

%% Question 4
%
% If you want to implement the non-linearity, you would have to keep track
% of the non-linear components using another vector. This would change the
% G matrix. This makes the time domain equation:
% 
% $$C\frac{dV}{dt}+GV+B=F$$
%
% Where B is the non-linear component-holding vector. In order to converge
% on a value for B, the simulation will have to not only iterate over time
% but also in each timestep.














