%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TI C2000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select numeric type, let's choose Q21
T=numerictype('WordLength',32,'FractionLength',21);
%Specify math attributes to the fimath object
F=fimath('RoundMode','floor','OverflowMode','wrap');
F.ProductMode='SpecifyPrecision';
F.ProductWordLength=32;
F.ProductFractionLength=21;
F.SumMode='SpecifyPrecision';
F.SumWordLength=32;
F.SumFractionLength=21;
%specify fipref object, to display warning in cases of overflow and
%underflow
P=fipref;
P.LoggingMode='on';
P.NumericTypeDisplay='none';
P.FimathDisplay='none';
%PLL Modelling starts from here
Fs=50000; %Sampling frequency = 50Khz
GridFreq=50; %Nominal Grid Frequency in Hz
Tfinal=0.2; %Time the simulation is run for = 0.5 seconds
Ts=1/Fs; %Sampling Time = 1/Fs
t=0:Ts:Tfinal; %Simulation Time vector
wn=2*pi*GridFreq; %Nominal Grid Frequency in radians
%generate input signal and create a fi object of it
%input wave with a phase jump at the mid point of simulation
% CASE 1 : Phase Jump at the Mid Point
L=length(t);
for n=1:floor(L)
u(n)=sin(2*pi*GridFreq*Ts*n);
end
for n=1:floor(L)
u1(n)=sin(2*pi*GridFreq*Ts*n);
end
for n=floor(L/2):L
u(n)=sin(2*pi*GridFreq*Ts*n+pi/2);
end
%CASE 2 : Harmonics
% L=length(t);
% for n=1:floor(L)
% u(n)=0.9*sin(2*pi*GridFreq*Ts*n)+0.1*sin(2*pi*5*GridFreq*Ts*n);
% end
% for n=1:floor(L)
% u1(n)=sin(2*pi*GridFreq*Ts*n);
% end
%CASE 3 : Frequency Shift
% L=length(t);
% for n=1:floor(L)
% u(n)=sin(2*pi*GridFreq*Ts*n);
% end
% for n=1:floor(L)
% u1(n)=sin(2*pi*GridFreq*Ts*n);
% end
% for n=floor(L/2):L
% u(n)=sin(2*pi*GridFreq*1.1*Ts*n);
% end
%CASE 4: Amplitude Variations
% L=length(t);
% for n=1:floor(L)
% u(n)=sin(2*pi*GridFreq*Ts*n);
% end
% for n=1:floor(L)
% u1(n)=sin(2*pi*GridFreq*Ts*n);
% end
% for n=floor(L/2):L
% u(n)=0.8*sin(2*pi*GridFreq*Ts*n);
% end;
u=fi(u,T,F);
u1=fi(u1,T,F);
%declare arrays used by the PLL process
Upd=fi([0,0,0],T,F);
ynotch=fi([0,0,0],T,F);
ynotch_buff=fi([0,0,0],T,F);
ylf=fi([0,0],T,F);
SinGen=fi([0,0],T,F);
Plot_Var=fi([0,0],T,F);
Mysin=fi([0,0],T,F);
Mycos=fi([fi(1.0,T,F),fi(1.0,T,F)],T,F);
theta=fi([0,0],T,F);
werror=fi([0,0],T,F);
%notch filter design
c1=0.1;
c2=0.00001;
X=2*c2*wn*2*Ts;
Y=2*c1*wn*2*Ts;
Z=wn*2*wn*2*Ts*Ts;
B_notch=[1 (X-2) (-X+Z+1)];
A_notch=[1 (Y-2) (-Y+Z+1)];
B_notch=fi(B_notch,T,F);
A_notch=fi(A_notch,T,F);
% simulate the PLL process
for n=2:Tfinal/Ts % No of iteration of the PLL process in the simulation time
% Phase Detect
Upd(1)= u(n)*Mycos(2);
%Notch Filter
ynotch(1)=-A_notch(2)*ynotch(2)-
A_notch(3)*ynotch(3)+B_notch(1)*Upd(1)+B_notch(2)*Upd(2)+B_notch(3)*Upd(3);
%update the Upd array for future sample
Upd(3)=Upd(2);
Upd(2)=Upd(1);
% PI Loop Filter
%ts=30ms, damping ration = 0.7
% we get natural frequency = 110, Kp=166.6 and Ki=27755.55
% B0=166.877556 & B1=-166.322444
ylf(1)= fi(1.0,T,F)*ylf(2)+fi(166.877556,T,F)*ynotch(1)+fi(-166.322444,T,F)*ynotch(2);
%update Ynotch for future use
ynotch(3)=ynotch(2);
ynotch(2)=ynotch(1);
ynotch_buff(n+1)=ynotch(1);
ylf(1)=min([ylf(1) fi(200.0,T,F)]);
ylf(2)=ylf(1);
wo=fi(wn,T,F)+ylf(1);
werror(n+1)=(wo-wn)*fi(0.00318309886,T,F);
%integration process
Mysin(1)=Mysin(2)+wo*fi(Ts,T,F)*(Mycos(2));
Mycos(1)=Mycos(2)-wo*fi(Ts,T,F)*(Mysin(2));
%limit the oscillator integrators
Mysin(1)=max([Mysin(1) fi(-1.0,T,F)]);
Mysin(1)=min([Mysin(1) fi(1.0,T,F)]);
Mycos(1)=max([Mycos(1) fi(-1.0,T,F)]);
Mycos(1)=min([Mycos(1) fi(1.0,T,F)]);
Mysin(2)=Mysin(1);
Mycos(2)=Mycos(1);
%update the output phase
theta(1)=theta(2)+wo*Ts;
%output phase reset condition
if(Mysin(1)>0 && Mysin(2) <=0)
theta(1)=-fi(pi,T,F);
end
SinGen(n+1)=Mycos(1);
Plot_Var(n+1)=Mysin(1);
end
% CASE 1 : Phase Jump at the Mid Point
error=Plot_Var-u;
%CASE 2 : Harmonics
%error=Plot_Var-u1;
%CASE 3: Frequency Variations
%error=Plot_Var-u;
%CASE 4: Amplitude Variations
%error=Plot_Var-u1;
figure;
subplot(3,1,1),plot(t,Plot_Var,'r',t,u,'b'),title('SPLL(red) & Ideal Grid(blue)');
subplot(3,1,2),plot(t,error,'r'),title('Error');
subplot(3,1,3),plot(t,u1,'r',t,Plot_Var,'b'),title('SPLL Out(Blue) & Ideal Grid(Red)');
