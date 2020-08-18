%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subject: Phase Synch. Evalulation
%          of Real Data (Kirby21)
% looking at the heatmap of time transition
% Author: Hamed Honari
% Advisor: M. Lindquist
% Date: NOV.10, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empirical Illustration of the Methods for Phase Synchronization
% Dataset: Kirby21

clear;clc;
load('K21_signchanged.mat');

N = 20;   % Number of subjects
nS = 2;   % number of states
nR = 2;   % number of runs
winLen = 30;   % Window length in Sliding Window

% 21 subjects 2 measurements (1 subject dropped cause it was bad)
ord{1} = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39];
ord{2} = [2 4 6 8 10 12 14 16 18 20 2 24 26 28 30 32 34 36 38 40];

% Removing the noisy ROIs from the data
for i = 1:1:nR*N
   Data{i} = Data2{i}(:,[32 25 7 13 23 2 4 5 14 9 12 16 19 20 26 29 3 10 27 31 24]); 
   %Data{i} = Data2{i}(:,[32 25 7 13 23 2 4 5 14 9 12 16 19 20 26 29 3 10 27 31 24 1 6 8  11 15 17 18 21 22 28 30 33 34 35 36 37 38 39]);
end

indx = nchoosek(1:21,2);
m = 1;

TR = 2;   % Kirby21 TR/TE = 2000/30 ms ==> TR = 2 [s]
Fs = 1/TR;


for m =1:N
signal = Data{1,m};
% Bandpass filtering the data
% the cutoff frequencies in bandpass should be in pi radians/sample
%filtered_signal = bandpass(signal,[0.03 0.07],Fs);

ftype = 'bandpass';

[z,p,kf] = butter(5,[0.03 0.07]/(Fs/2),ftype);
[b,a] = zp2tf(z,p,kf);
filtered_signal=filtfilt(b,a,signal);
hilb_filtered_signal = hilbert(filtered_signal);
phi{m} = angle(hilb_filtered_signal);

for j = 1:size(Data{1,1},1)
    for i = 1:1:size(indx,1)
        delphi{m}(j,i) = bsxfun(@minus, phi{m}(j,indx(i,1)),phi{m}(j,indx(i,2)));
    end
end

for i = 1:size(indx,1)
PLV{m,1}(:,i) = phaseloc(phi{m}(:,indx(i,1)),phi{m}(:,indx(i,2)),'option','window','winsize',winLen);
CCORSW{m,1}(:,i) = circcorSW(phi{m}(:,indx(i,1)),phi{m}(:,indx(i,2)),'option','window','winsize',winLen);
NCCORSW{m,1}(:,i) = newcirccorSW(phi{m}(:,indx(i,1)),phi{m}(:,indx(i,2)),'option','window','winsize',winLen);
SINDELPHI1{m}(:,i)= 1-abs(sin(phi{m}(:,indx(i,1))-phi{m}(:,indx(i,2))));
COSDELPHI1{m}(:,i) = cos(phi{m}(:,indx(i,1))-phi{m}(:,indx(i,2)));
end

smoothSin{m} = movmean(SINDELPHI1{m},winLen,1,'Endpoints','fill');
smoothCos{m} = movmean(COSDELPHI1{m},winLen,1,'Endpoints','fill');
m
end
%% ---------------------------------------------------


display('Dont worry about the warning, it just ignores the ends in sliding window')


% k-means clustering of the matrices of phase synch.
[idx{1},Corr{1}] = mykmeans(cat(1,PLV{:,1}),nS); 
[idx{2},Corr{2}] = mykmeans(cat(1,CCORSW{:,1}),nS);
[idx{3},Corr{3}] = mykmeans(cat(1,NCCORSW{:,1}),nS);
[idx{4},Corr{4}] = mykmeans(cat(1,SINDELPHI1{:}),nS);
[idx{5},Corr{5}] = mykmeans(cat(1,COSDELPHI1{:}),nS);
% [idxsmoothsin,Corrsmoothsin] = mykmeans(cat(1,smoothSin{:}),nS);
% [idxsmoothcos,Corrsmoothcos] = mykmeans(cat(1,smoothCos{:}),nS);


% Using the sliding window
SWDat = cat(1,Data{1,1:20});
[CSW,~] = Sliding_Window(SWDat,winLen);
for i=1:size(CSW,3)
    CSWVect(i,:) = vmconv(CSW(:,:,i),'mat2vec');
end
[idx{6},Corr{6}] = mykmeans(CSWVect,nS);
% Sliding Window Pre-Whitening using AR(1)
 [CV2] = sliding_window_white_AR1(SWDat, winLen,'white');
for i=1:size(CV2,3)
    CV2Vect(i,:) = vmconv(CV2(:,:,i),'mat2vec');
end
[idx{7},Corr{7}] = mykmeans(CV2Vect,nS);




count = 0;
figure;
h = tight_subplot(7, nS, [.001 .001],[.01 .001],[.05 .05]);
labels = {'PLV','$\rho_{circ}$','$\rho_{tor}$','$\Psi$','$\vartheta$','CSW','PW-CSW'};
for i = 1:7
    matchix = matchstates(Corr{1},Corr{i},1); % to match the states
    for j = 1:nS
        count = count + 1;
        axes(h(count));gsplot(Corr{i}(:,:,matchix(j)));
        axis square;
        set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w')
        if j == 1
            ylabel(labels{i},'interpreter','latex');
        end
    end
    
end

suptitle('Kirby21-Dataset')
