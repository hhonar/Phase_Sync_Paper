%%%%%%%%%%%%%%%%%%%%%
% Script for Simulations of Phase Synchronization [Simulation 2 & 3]
% Written by: Hamed Honari
% Advisor: Martin A. Lindquist
%%%%%%%%%%%%%%%%%%%%%
% Dependencies (used subroutines: CircularStat.m, phaseloc.m, circcorSW,...
% ncirccorSW, boundedline, ipsprep)

clear;clc
% Let's assume the situation where TR is 2 [s] similar to the Kirby21 data
TR = 2;                                                 % Repetition time 
fs = 1/TR;                                           % Sampling frequency
t = 0:1/fs:668-1/fs;
% Defining freq. component of signals to be generated for the simulations
f = 0.05; % freq. component of the signal x1

smltn = input('Which simulation do you want to run: 2. Ramp 3. Sigmoid?')
switch smltn
    case 2
        phi2 = 4*pi/334.*(t-334).*(t-334>=0);
    case 3
        phi2 = 2*pi./(1+exp(-0.01*(t-334)));
end

% Creating the sinusoidal signals
x1 = cos(2*pi*f*t);        % first signal
y1 = cos(2*pi*f*t + phi2); % second signal


% Plotting the signals
figure;hold on;
plot(t,x1,'r');
plot(t,y1,'b','LineWidth',2);
xlabel('time');ylabel('Amplitude [au]');
title('Simulated sinusoidal signals without additive white noise')

figure;periodogram(x1,[],'onesided',2^nextpow2(length(x1)),fs)
hold on;periodogram(y1,[],'onesided',2^nextpow2(length(y1)),fs)

% Plotting the signals
figure;hold on;
plot(t,x1,'r','LineWidth',2);
plot(t,y1,'b','LineWidth',2);
xlabel('time');ylabel('Amplitude [au]');
%title('Simulated sinusoidal signals without additive white noise')
xticklabels(1:334)

method = input('Do you want to bandpass filter? type 1 for no, 2 for yes');
switch method
    case 1
        strg = 'nofilter';
    case 2
        strg = 'filtfiltbutw';
end
w = [30,60,120];

for m = 1:1:200
noise = mvnrnd([0 0],[1 0;0 1],length(t))';
ex = noise(1,:);
ey = noise(2,:);
XN1 = x1 + ex;
YN1 = y1 + ey;
Data = [XN1;YN1];
% calculating the PLV for the pair of signals (using the written m-script
% 'phaseloc.m'
[IPData1] = ipsprep(Data,0.03,0.07,strg,fs,5,3,100);
DELPHI{m} = IPData1(:,1)-IPData1(:,2);
PLV{1}(:,m) = phaseloc(IPData1(:,1),IPData1(:,2),'option','window','winsize',w(1));
PLV{2}(:,m) = phaseloc(IPData1(:,1),IPData1(:,2),'option','window','winsize',w(2));
PLV{3}(:,m) = phaseloc(IPData1(:,1),IPData1(:,2),'option','window','winsize',w(3));
CCORSW{1}(:,m) = circcorSW(IPData1(:,1),IPData1(:,2),'option','window','winsize',w(1));
CCORSW{2}(:,m) = circcorSW(IPData1(:,1),IPData1(:,2),'option','window','winsize',w(2));
CCORSW{3}(:,m) = circcorSW(IPData1(:,1),IPData1(:,2),'option','window','winsize',w(3));
NCCORSW{1}(:,m) = newcirccorSW(IPData1(:,1),IPData1(:,2),'option','window','winsize',w(1));
NCCORSW{2}(:,m) = newcirccorSW(IPData1(:,1),IPData1(:,2),'option','window','winsize',w(2));
NCCORSW{3}(:,m) = newcirccorSW(IPData1(:,1),IPData1(:,2),'option','window','winsize',w(3));
SINDELPHI1(:,m)= 1-abs(sin(IPData1(:,1)-IPData1(:,2)));
COSDELPHI1(:,m) = cos(IPData1(:,1)-IPData1(:,2));
m
end



%% Displaying the results

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure;
subplot(6,1,1);hold on;plot((1:length(t)),phi2,'k','LineWidth',1.5);box on;
xlabel('time [s]','interpreter','latex');
ylabel('$\Delta\Phi[t] \quad [rad]$','interpreter','latex');
title('(a)','interpreter','latex');
legend('$\Delta\Phi[t]$','Location','northwest');xlim([0 450]);
subplot(6,1,2);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(PLV{1},2),0.95.*std(PLV{1},1,2), '-r','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
subplot(6,1,2);hold on;[hl2 hp2]=boundedline((1:length(t)),mean(PLV{2},2),0.95.*std(PLV{2},1,2), '-b','alpha','nan','remove');
outlinebounds(hl2,hp2)
subplot(6,1,2);hold on;[hl3 hp3]=boundedline((1:length(t)),mean(PLV{3},2),0.95.*std(PLV{3},1,2), '-g','alpha','nan','remove');box on;
outlinebounds(hl3,hp3)
legend([hp1 hp2 hp3],strcat(['Window length  = ' num2str(w(1))]),strcat(['Window length  = ' num2str(w(2))]),strcat(['Window length  = ' num2str(w(3))]),'Location','southeast');
xlabel('time [s]','interpreter','latex');
ylabel('$PLV$','interpreter','latex');
title('(c)','interpreter','latex');
set([hl1 hl2 hl3],'LineWidth',2);xlim([0 450]);ylim([0 1])
subplot(6,1,3);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(CCORSW{1},2),0.95.*std(CCORSW{1},1,2), '-r','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
subplot(6,1,3);hold on;[hl2 hp2]=boundedline((1:length(t)),mean(CCORSW{2},2),0.95.*std(CCORSW{2},1,2), '-b','alpha','nan','remove');box on;
outlinebounds(hl2,hp2)
subplot(6,1,3);hold on;[hl3 hp3]=boundedline((1:length(t)),mean(CCORSW{3},2),0.95.*std(CCORSW{3},1,2), '-g','alpha','nan','remove');box on;
outlinebounds(hl3,hp3)
set([hl1 hl2 hl3],'LineWidth',2)
legend([hp1 hp2 hp3],strcat(['Window length  = ' num2str(w(1))]),strcat(['Window length  = ' num2str(w(2))]),strcat(['Window length  = ' num2str(w(3))]),'Location','southeast');
xlabel('time [s]','interpreter','latex');
ylabel('$\rho_{circ}$','interpreter','latex');
title('(d)','interpreter','latex');xlim([0 450]);ylim([-1 1])
subplot(6,1,4);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(NCCORSW{1},2),0.95.*std(NCCORSW{1},1,2), '-r','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
subplot(6,1,4);hold on;[hl2 hp2]=boundedline((1:length(t)),mean(NCCORSW{2},2),0.95.*std(NCCORSW{2},1,2), '-b','alpha','nan','remove');box on;
outlinebounds(hl2,hp2)
subplot(6,1,4);hold on;[hl3 hp3]=boundedline((1:length(t)),mean(NCCORSW{3},2),0.95.*std(NCCORSW{3},1,2), '-g','alpha','nan','remove');box on;
outlinebounds(hl3,hp3)
set([hl1 hl2 hl3],'LineWidth',2)
legend([hp1 hp2 hp3],strcat(['Window length  = ' num2str(w(1))]),strcat(['Window length  = ' num2str(w(2))]),strcat(['Window length  = ' num2str(w(3))]),'Location','southeast');
xlabel('time [s]','interpreter','latex');
ylabel('$\rho_{tor}$','interpreter','latex');
title('(e)','interpreter','latex');xlim([0 450]);ylim([-1 1])
subplot(6,1,5);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(SINDELPHI1,2),0.95.*std(SINDELPHI1,1,2), '-r','alpha','nan','remove');box on;
legend('$1-|sin(\Delta\Phi[t])|$','Location','Best');
xlabel('time [s]','interpreter','latex');
ylabel('$\Psi$','interpreter','latex');
title('(f)','interpreter','latex');xlim([0 450]);ylim([0 1])
subplot(6,1,6);hold on;[hl1 hp1]=boundedline((1:length(t)),mean(COSDELPHI1,2),0.95.*std(COSDELPHI1,1,2), '-m','alpha','nan','remove');box on;
legend('$cos(\Delta\Phi[t])$','Location','Best');
xlabel('time [s]','interpreter','latex');
ylabel('$\vartheta$','interpreter','latex');
title('(g)','interpreter','latex');xlim([0 450]);ylim([-1 1])




