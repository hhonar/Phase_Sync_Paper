%%%%%%%%%%%%%%%%%%%%%
% Script for Bivariate Signals - Null Setting WITH SURROGATE
% Written by: Hamed Honari
% Advisor: Martin A. Lindquist
%%%%%%%%%%%%%%%%%%%%%
% Dependencies (used subroutines: CircularStat.m, phaseloc.m, circcorSW,...
% ncirccorSW, boundedline, ipsprep)

clear;clc;
% Let's assume the situation where TR is 2 [s] similar to the Kirby21 data
TR = 2; % Repetition time 
fs = 1/TR; % Sampling frequency
t = 0:1/fs:1000-1/fs;
% Defining the freq. component of the signal
f = 0.05; % freq. component of the signal x1



Data1 = mvnrnd([0 0],[1 0;0 1],334)';
x1 = Data1(1,:);
y1 = Data1(2,:);

%y1 = wgn(1,length(t),0);
%[surrX,params]=surrogate(x1, 1000, 'TS', 0, fs);

% Option 1: commenting it out 
% for m = 1:1:1000    
% surrX(m,:) = osl_surrogate(x1','phase_randomization')';
% surrY(m,:) = osl_surrogate(y1','phase_randomization')';
% end


N = 1000;

% Creating the surrogates
[surrX,paramsX] = surrogate(x1, N, 'RP', 1, fs);
[surrY,paramsY] = surrogate(y1, N, 'RP', 1, fs);

% Ensuring the length of the surrogates created for signals are the same
while (size(surrX,2) ~= size(surrY,2))
    Data1 = mvnrnd([0 0],[1 0;0 1],334)';
    x1 = Data1(1,:);
    y1 = Data1(2,:);
    [surrX,paramsX] = surrogate(x1, 1000, 'RP', 1, fs);
    [surrY,paramsY] = surrogate(y1, 1000, 'RP', 1, fs);
end



method = input('Do you want to bandpass filter? type 1 for no, 2 for yes');
switch method
    case 1
        strg = 'nofilter';
    case 2
        strg = 'filtfiltbutw';
end

w = [30,60,120];

for m = 1:1:N
X1 = surrX(m,:); 
Y1 = surrY(m,:); 
Data = [X1;Y1];
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
end





%% Display the results

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure;
subplot(5,1,1);hold on;[hl2 hp2]=boundedline((1:length(IPData1)),mean(PLV{1},2),0.95.*std(PLV{1}'), '-r','alpha','nan','remove');box on;
outlinebounds(hl2,hp2)
subplot(5,1,1);hold on;[hl3 hp3]=boundedline((1:length(IPData1)),mean(PLV{2},2),0.95.*std(PLV{2}'), '-b','alpha','nan','remove');
outlinebounds(hl3,hp3)
subplot(5,1,1);hold on;[hl4 hp4]=boundedline((1:length(IPData1)),mean(PLV{3},2),0.95.*std(PLV{3}'), '-g','alpha','nan','remove');box on;
outlinebounds(hl4,hp4)
legend([hp2 hp3 hp4],strcat(['Window length  = ' num2str(w(1))]),strcat(['Window length  = ' num2str(w(2))]),strcat(['Window length  = ' num2str(w(3))]),'Location','southeast');
xlabel('t [s]','interpreter','latex');
ylabel('$PLV$','interpreter','latex');
title('(a)','interpreter','latex');
set([hl4 hl2 hl3],'LineWidth',2);xlim([0 400]);ylim([0 1]);
subplot(5,1,2);hold on;[hl1 hp1]=boundedline((1:length(IPData1)),mean(CCORSW{1},2),0.95.*std(CCORSW{1},1,2), '-r','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
subplot(5,1,2);hold on;[hl2 hp2]=boundedline((1:length(IPData1)),mean(CCORSW{2},2),0.95.*std(CCORSW{2},1,2), '-b','alpha','nan','remove');box on;
outlinebounds(hl2,hp2)
subplot(5,1,2);hold on;[hl3 hp3]=boundedline((1:length(IPData1)),mean(CCORSW{3},2),0.95.*std(CCORSW{3},1,2), '-g','alpha','nan','remove');box on;
outlinebounds(hl3,hp3)
set([hl1 hl2 hl3],'LineWidth',2)
legend([hp1 hp2 hp3],strcat(['Window length  = ' num2str(w(1))]),strcat(['Window length  = ' num2str(w(2))]),strcat(['Window length  = ' num2str(w(3))]),'Location','southeast');
xlabel('t [s]','interpreter','latex');
ylabel('$\rho_{circ}$','interpreter','latex');
title('(b)','interpreter','latex');xlim([0 400]);ylim([-1 1]);
subplot(5,1,3);hold on;[hl1 hp1]=boundedline((1:length(IPData1)),mean(NCCORSW{1},2),0.95.*std(NCCORSW{1},1,2), '-r','alpha','nan','remove');box on;
outlinebounds(hl1,hp1)
subplot(5,1,3);hold on;[hl2 hp2]=boundedline((1:length(IPData1)),mean(NCCORSW{2},2),0.95.*std(NCCORSW{2},1,2), '-b','alpha','nan','remove');box on;
outlinebounds(hl2,hp2)
subplot(5,1,3);hold on;[hl3 hp3]=boundedline((1:length(IPData1)),mean(NCCORSW{3},2),0.95.*std(NCCORSW{3},1,2), '-g','alpha','nan','remove');box on;
outlinebounds(hl3,hp3)
set([hl1 hl2 hl3],'LineWidth',2)
legend([hp1 hp2 hp3],strcat(['Window length  = ' num2str(w(1))]),strcat(['Window length  = ' num2str(w(2))]),strcat(['Window length  = ' num2str(w(3))]),'Location','southeast');
xlabel('t [s]','interpreter','latex');
ylabel('$\rho_{tor}$','interpreter','latex');
title('(c)','interpreter','latex');xlim([0 400]);ylim([-1 1]);
subplot(5,1,4);hold on;[hl1 hp1]=boundedline((1:length(IPData1)),mean(SINDELPHI1,2),0.95.*std(SINDELPHI1,1,2), '-r','alpha','nan','remove');box on;
legend('$1-|sin(\Delta\Phi(t))|$','Location','Best');
xlabel('t [s]','interpreter','latex');
ylabel('$\Psi$','interpreter','latex');
title('(d)','interpreter','latex');xlim([0 400]);ylim([0 1]);
set([hl1],'LineWidth',2);
subplot(5,1,5);hold on;[hl1 hp1]=boundedline((1:length(IPData1)),mean(COSDELPHI1,2),0.95.*std(COSDELPHI1,1,2), '-m','alpha','nan','remove');box on;
legend('$cos(\Delta\Phi(t))$','Location','Best');
xlabel('t [s]','interpreter','latex');
ylabel('$\vartheta$','interpreter','latex');
title('(e)','interpreter','latex');xlim([0 400]);ylim([-1 1]);
set([hl1],'LineWidth',2)
set(findall(gcf,'-property','FontSize'),'FontSize',12)


