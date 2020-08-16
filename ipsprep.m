function [IPData1,filtsig1] = ipsprep(fMRI_data,fcutlow,fcuthigh,method,fs,order,Rp,Rs)



% Defining the range for the high cut-off frequency
%fcuthigh = 0.055;


DELT = length(fMRI_data);%2240;%:-10:250;

%method = 'connfiltcanlab';
%method = 'filtfiltbutw';


for j = 1:1:length(DELT)
for k = 1:1:1   % k index for changing the various high cut-off freq.

        %DATA1 = [X1(1,1:DELT(j));Y1(1,1:DELT(j))];
               
        %fs=1/30; %sample rate in Hz
        inputsig1=fMRI_data';  %input signal
        %order=4;   %order of filter (controls the ripples)
       % fcutlow=0.045;   %low cut frequency in Hz
    switch lower(method)
        case 'matbandpass'
            filtsig1 = bandpass(inputsig1,[fcutlow fcuthigh],fs,'ImpulseResponse','iir','Steepness',0.95);%,'ImpulseResponse','iir','Steepness',0.95
        case 'spmbandpass'
            % Using the SPM_FILTER function
K{1}.RT     = 1/fs;                  % assign TR to K
K{1}.HParam = 1./fcutlow;                   % assign HPF cutoff (s) to K
K{1}.row    = 1:length(inputsig1');           % the timepoints to filter (all, in this case)
K{1}.LParam = 1/fcuthigh;  
K{1}.HChoice = 'specify';
K{1}.LChoice = 'hrf';
K = pr_spm_filter('set'  ,K);         % add the filter matrix to the sturcture K  
fx = pr_spm_filter('high' ,K,inputsig1);        % filter the timeseries using K

% Performing the low-pass filter
clear K;
K{1}.RT     = 1/fs;                  % assign TR to K
K{1}.HParam = 1/fcutlow;                   % assign HPF cutoff (s) to K
K{1}.row    = 1:length(inputsig1');           % the timepoints to filter (all, in this case)
K{1}.LParam = 0.5/fcuthigh;  
K{1}.HChoice = 'specify';
K{1}.LChoice = 'Gaussian';
K = pr_spm_filter('set'  ,K);         % add the filter matrix to the sturcture K  
ffx = pr_spm_filter('low' ,K,fx);        % filter the timeseries using K
filtsig1 = ffx;            
            
            
        case 'filtfiltbutw'
        %fcuthigh=0.06;   %high cut frequency in Hz
            ftype = 'bandpass';
            %[b,a]=butter(order,[fcutlow fcuthigh]/(fs/2));
            %filtsig1=filtfilt(b,a,inputsig1);  %filtered signal
            
            [z,p,kf] = butter(order,[fcutlow fcuthigh]/(fs/2),ftype);
            %sos = zp2sos(z,p,k);
            [b,a] = zp2tf(z,p,kf);
            filtsig1=filtfilt(b,a,inputsig1);
        case 'filtfiltelip' % NOTE: if order is assigned 10 indeed the actual order is 20
            ftype = 'bandpass';
            [z,p,kf] = ellip(order,Rp,Rs,[fcutlow fcuthigh]/(fs/2),ftype);
            [b,a] = zp2tf(z,p,kf);
            filtsig1=filtfilt(b,a,inputsig1);
        case 'connfiltcanlab'
            filtsig1 = wgr_band_filter(inputsig1,1/fs,[fcutlow fcuthigh]);
        case 'serj'
            filtsig1 = BPFilter5(inputsig1,fs,[fcutlow fcuthigh],order);
        case 'nofilter'
            filtsig1 = inputsig1;
    end
        %[filtsig1, ~, ~, ~] = hpfilter(y, 1/fs, [0.045 0.055], spersess);
% Extracting the instantenous phase by calling the subroutine 
        [IPData1, IP_matrixData1,~,dPhiData1{j}] = iphsync(filtsig1);
        
        j

end

DELTAPHI = IPData1(:,1) - IPData1(:,2);
clear i;
RX1 = exp(i*DELTAPHI(:,1));

end


end %of main function

%% ------------------------------------------------------------------------
% S U B R O U T I N E
%% ------------------------------------------------------------------------

% Instantaneous Phase Synchrony Function
function [IP, IP_matrix, IP_matrix1, dPhiXY] = iphsync(fMRI_data)
fprintf('\n \t Calculating Instantaneous Phase Synchrony (IPS) ...\n')
reverseStr = '';
theElapsedTime = tic;

%[~, IP] = hilbertFIR(fMRI_data); 

Z(:,1) = hilbert(fMRI_data(:,1));
Z(:,2) = hilbert(fMRI_data(:,2));
IP = angle(Z);

%IP = angle(hilbert(fMRI_data)); % instantaneous phase of NxT fMRI data; node(N) by time (T).

IP_matrix = zeros(size(IP,1), size(IP,1), size(IP,2)); % Pre-allocate 3D matrix

for time_point = 1:size(IP,2) % Total number of fMRI time-points
    IP_matrix(:,:,time_point) = 1 - abs(sin(bsxfun(@minus,IP(:,time_point)', IP(:,time_point)))); % node-by-node synchrony at each time point
    IP_matrix1(:,:,time_point)= bsxfun(@minus,IP(:,time_point)', IP(:,time_point));
    msg = sprintf('\n \t IPS time-point %d/%d ...', time_point, size(IP,2)); % on screen information of calculation progress
    fprintf([reverseStr,msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
end

dPhiXY = 1-abs(sin(IP(:,1)-IP(:,2)));
theElapsedTime = toc(theElapsedTime);
fprintf('\n\t Elapsed time: %g seconds ...\n', theElapsedTime);

%figure;plot(fMRI_data');
%figure;histogram(IP_matrix(:,:,1),'FaceAlpha',0.2)
%hold on;histogram(IP_matrix(:,:,2),'FaceColor','r','FaceAlpha',0.3);
%xlabel('')
%figure;histogram(IP(:,1),40,'FaceAlpha',0.7)
%hold on;histogram(IP(:,2),40,'FaceColor','r','FaceAlpha',0.3)

end
