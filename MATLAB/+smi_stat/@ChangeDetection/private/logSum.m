% Mark J. Olah [mjo@cs.unm.edu]
% 2014
function logS = logSum(varargin)
    % This allows for efficient computation of the log of the sum
    % of a bunch of values that are only known in logarithm.
    % Returns logS = log(sum(exp(logX)))
    if nargin==2
        log_a = max(varargin{1},varargin{2});
        delta = abs(varargin{1}-varargin{2});
        if delta > 20
            logS = log_a + log1p(exp(-delta));
        else
            logS = log_a + log(1+exp(-delta));
        end        
    else
        if nargin==1
            logX=varargin{1};
        else
            logX=[varargin{:}];
        end
        C = max(logX);
        logX = logX-C;
        MIN_LOG = -1e2;
        logS = C + log(sum(exp( logX(logX>MIN_LOG) )));
    end
end
