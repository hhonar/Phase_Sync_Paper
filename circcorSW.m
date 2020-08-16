%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for Sliding Window using Circular-Circular Correlation
% Sliding window of the Phase Locking Value
% Date: Dec. 3, 2018
% Author: Hamed Honari
% Advisor: Prof. Martin A. Lindquist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [rho] = circcorSW(alpha, beta, varargin)
% This function requires 2 inputs and has some optional inputs.
% Inputs:
%         alpha: (extracted) phase-time series of signal 1
%         beta:  (extracted) phase-time series of signal 2 
% only want 2 optional inputs at most
%         opt1: could be 'whole' (default),'window'. 
%         opt2: could be 'boxcar' (default),'hamming'.

% First Checking the dependent script - CircStatToolbox
toolCheck = exist('circstat-matlab-master');

switch toolCheck
    case 7
        % DO NOTHING! The toolbox is already installed!
    case 0 
        % INSTALL the circ-stat toolbox developed by Philipp Berens from
        % Github
        % Getting the search path of the MATLAB
        file = path;
        [filepath,~,~] = fileparts(file);
        searchpath = regexp(filepath,':','split'); 
        searchpath{1};
        currentFolder = pwd;
        cd(searchpath{1});
        disp('NOTE: This will install the circstat toolbox for Circular Statistics in your search path since it is a dependency');
        url = 'https://github.com/circstat/circstat-matlab/archive/master.zip';
        unzip(url,pwd);
        cd(currentFolder);
end

p = inputParser;
addRequired(p,'alpha',@isvector);
addRequired(p,'beta',@isvector);

defaultWinSize = 0;   % add the optional window size  with this default value
addOptional(p,'winsize',defaultWinSize,@isnumeric)

defaultOption = 'whole';
checkString = @(s) any(strcmp(s,{'whole','window'}));
addParameter(p,'option',defaultOption,checkString);

defaultWinType = 'boxcar';
checkString1 = @(s) any(strcmp(s,{'boxcar','hamming'}));
addParameter(p,'wintype',defaultWinType,checkString1);

parse(p,alpha,beta,varargin{:});
alpha   = p.Results.alpha;
beta    = p.Results.beta;
winsize = p.Results.winsize;
option  = p.Results.option;
wintype = p.Results.wintype;
%numvarargs = nargin;
%if numvarargs > 3
%    error('Error in phaseloc:TooManyInputs', ...
%        'requires at most 3 optional inputs');
%end

% set defaults for optional inputs
%optargs = {'whole' 'boxcar' 10};
%optargs = {1 1 10};


% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
%optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
%[opt1, opt2, win] = optargs{:};
rho = NaN(size(alpha)); % allocate 

switch option
    case {'whole'}
        clear circCorr;
        rho = circ_corrcc(alpha,beta);
    case {'window'}
        if winsize ~= 0
            switch wintype
                case {'boxcar'}
                    for k= winsize:length(alpha)
                        rho(k-winsize/2) = circ_corrcc(alpha(k-winsize+2:k-1),beta(k-winsize+2:k-1));
                    end   
                case {'hamming'}
                    for k= winsize:length(alpha)
                        rho(k-winsize/2) = circ_corrcchamming(alpha(k-winsize+1:k),beta(k-winsize+1:k));
                    end
            end
        elseif winsize == 0
            clear circCorr;
            rho = circ_corrcc(alpha,beta);
            disp('Warning: window size of 0 was assigned which is equavalently interpreted as no windowing option [static].  Hence the analysis was done on whole')
        end
end
        

end