%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for Sliding Window using a newly defined Circular-
% Circular Correlation (see Zhan et al. 2017)
% Date: Dec. 3, 2018
% Author: Hamed Honari
% Advisor: Prof. Martin A. Lindquist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [rho] = newcirccorSW(alpha, beta, varargin)
% This function requires 2 inputs and has some optional inputs.
% Inputs:
%         alpha: (extracted) phase-time series of signal 1
%         beta:  (extracted) phase-time series of signal 2 
% only want 2 optional inputs at most
%         opt1: could be 'whole' (default),'window'. 
%         opt2: could be 'boxcar' (default),'hamming'.

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
        rho = newcirccorr(alpha,beta);
    case {'window'}
        if winsize ~= 0
            switch wintype
                case {'boxcar'}
                    for k= winsize:length(alpha)
                        rho(k-winsize/2) = newcirccorr(alpha(k-winsize+1:k),beta(k-winsize+1:k));
                    end   
                case {'hamming'}
                    for k= winsize:length(alpha)
                        rho(k-winsize/2) = newcirccorrhamming(alpha(k-winsize+1:k),beta(k-winsize+1:k));
                    end 
            end
        elseif winsize == 0
            clear circCorr;
            rho = newcirccorr(alpha,beta);
            disp('Warning: window size of 0 was assigned which is equavalently interpreted as no windowing option.  Hence the analysis was done on whole')
        end
end
        

end


%------------------------------------------------
% SUBROUTINE 1: Order function h(alpha,beta)
%------------------------------------------------

% According to the paper the oder function is defined as:
% h(alpha,beta) = [(alpha-beta+2pi) mod 2pi] - pi
function [h] = orderfnc(alpha,beta)
h = mod((alpha - beta + 2*pi),2*pi) - pi;

%% Sanity check for h: PASSED 
% delta = alpha - beta;
% for i = 1:length(delta)
% if (delta(i) > -2*pi && delta(i) < 0)
%     f(i) = delta(i) + pi;
% elseif (delta(i) >= 0 && delta(i) < 2*pi)
%     f(i) = delta(i) - pi;
% end


end



%---------------------------------------------------
% SUBROUTINE 2: Estimation of new circ. correlation
%---------------------------------------------------

function rhohat = newcirccorr(alpha,beta)
    v = uint16(1:length(alpha));
    C = nchoosek(v,uint16(2));
    rhonumer = sum( orderfnc( alpha(C(:,1)),alpha(C(:,2)) ) .* orderfnc( beta(C(:,1)),beta(C(:,2)) )); 
    rhodenom = sqrt(sum((orderfnc(alpha(C(:,1)),alpha(C(:,2)))).^2).*sum((orderfnc(beta(C(:,1)),beta(C(:,2)))).^2));
    rhohat = rhonumer/rhodenom;
end



%---------------------------------------------------
% SUBROUTINE 3: Estimation of new circ. correlation - tapered
%---------------------------------------------------

function rhohat = newcirccorrhamming(alpha,beta)
    v = uint16(1:length(alpha));
    C = nchoosek(v,uint16(2));
    rhonumer = sum( hamming(length(alpha)).*(orderfnc( alpha(C(:,1)),alpha(C(:,2)) ) .* orderfnc( beta(C(:,1)),beta(C(:,2)) ))); 
    rhodenom = sqrt(sum(hamming(length(alpha)).*((orderfnc(alpha(C(:,1)),alpha(C(:,2)))).^2)).*sum(hamming(length(alpha)).*((orderfnc(beta(C(:,1)),beta(C(:,2)))).^2)));
    rhohat = rhonumer/rhodenom;
end


function rhohat = newcirccorr1(alpha,beta)
    v = uint16(1:length(alpha)-1);
    C = [v' (v+1)'];
    rhonumer = sum( orderfnc( alpha(C(:,1)),alpha(C(:,2)) ) .* orderfnc( beta(C(:,1)),beta(C(:,2)) )); 
    rhodenom = sqrt(sum((orderfnc(alpha(C(:,1)),alpha(C(:,2)))).^2).*sum((orderfnc(beta(C(:,1)),beta(C(:,2)))).^2));
    rhohat = rhonumer/rhodenom;
end
