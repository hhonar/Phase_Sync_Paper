function [ Ct,whitenedData ] = Sliding_Window(dat,windowsize,varargin)
% function [ output ] = sliding_window(,windowsize)
%
% Sliding window correlation of data
% 
% INPUTS:
%
%      dat          Zero mean T by p matrix
% 
% OUTPUTS:
%
%      Ct           p by p by T array of conditional correlations
%
%
% File created by Martin Lindquist on 07/22/14
% Last update: 07/22/14

[T,p] = size(dat);

% initalize optional variables to default values here.

dowhiten = 0;

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'white', 'whiten'}, dowhiten = 1;
                
            case {'whiteARParam'}, dowhiten = 2;
                                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if dowhiten == 1 
    
    Mdl = arima(1,0,0);
    
    opts = optimoptions('fmincon');
    opts.Algorithm = 'sqp';
    deltaX = 10^(-5);
    opts.TolX = deltaX;
    opts.TolCon = 10^(-6);
    opts.MaxFunEvals = 1000;
    
    indx = 0;
   % fprintf('ARMA(1,1) on %d variables: %05d', p, indx);

    for k=1:p
        
        EstMdl = estimate(Mdl, dat(:,k), 'Display', 'off', 'options', opts);
        
        [res, v] = infer(EstMdl,dat(:,k));
        dat(:,k) = res./sqrt(v);
        
        indx = indx + 1;
        fprintf('\b\b\b\b\b%05d', indx);
    end
    
    fprintf('\n');
    
elseif dowhiten == 2
    % assigning the AR(1) parameter \phi
    f(1) = 1;f(2) = -varargin{i};
    disp('whitening with specificied \phi parameter')
    for k=1:1:p
        dat(:,k) = filter(f,1,dat(:,k));
    end
end




Ct = NaN(p,p,T);
whitenedData = dat;

for i=(windowsize):T
    Ct(:,:,i)=corr(dat((i-windowsize+1):i,:));
end

end
