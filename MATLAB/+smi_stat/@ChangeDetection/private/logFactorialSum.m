% Mark J. Olah (mjo@cs.unm.edu)
%
% 06-2014
%
function LF=logFactorialSum(data)
    %
    % Compute the sum of the log of the factorial of all value in an array of values.
    % For numbers less than the threshold, we compute the factorial directly
    % For larger values we use stirlings approximation with the first order
    % correction term.
    %
    % Inputs:
    %   data: vector of integers to take log factorial of
    %
    % Outputs:
    %   LF: A scalar = sum(log(factorial(data))).
    threshold=30;
    bigvals=data(data>threshold);
    smallvals=data(data<=threshold);
    LFsmall=sum(log(factorial(smallvals)));
    LFbig=sum(bigvals.*log(bigvals)-bigvals+0.5*log(2*pi*bigvals));
    LF=LFbig+LFsmall;
end
