% Mark J. Olah (mjo@cs.unm.edu)
%
% 06-2014
%

function LF=logFactorial(data)
    %
    % Compute the log of the factorial of each value in an array of values.
    % For numbers less than the threshold, we compute the factorial directly
    % For larger values we use stirlings approximation with the first order
    % correction term.
    %
    % Inputs:
    %   data: vector of integers to take log factorial of
    %
    % Outputs:
    %   LF: A vector of the same size as the inputs where each value is
    %       log(data(i)!).
    threshold=30;
    bigvals=data(data>threshold);
    smallvals=data(data<=threshold);
    LFsmall=log(factorial(smallvals));
    LFbig=bigvals.*log(bigvals)-bigvals+0.5*log(2*pi*bigvals);
    LF=zeros(size(data));
    LF(data<=threshold)=LFsmall;
    LF(data>threshold)=LFbig;
end
