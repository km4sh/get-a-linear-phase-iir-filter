%% Equiripple Group Delay Error Allpass Filters Designing
%  log:
%   version1:
%    july 10, 2017
%   version2:(nightly)
%    july 11, 2017
%      1. added the non-linear constraint.
%      2. able to work out some coefficient.
%
%  issues:
%   options trying:
%    StepTolerance;
%    FunctionTolerance;
%
%
%  author:
%   matt ma @scie, shanghai university
%   mattma9209@aliyun.com
%
%  input vars:
%   freqp: frequency points;
%   phred: phase response desired.
%
%  output vars:
%   coeff: coefficient vector of designed allpass filter as desired.
%
%  reference:
%   Design and Application of Allpass Filters with Equiripple Group Delay Errors(2013);

function coeff = eqrpgdr(freqp, phred, order)

    %% step 1 initialize
    coeff = zeros(1,order);
    olderr = inf;
    delta = 1e-8;
    coeffinit = zeros(1,order); % correct?
    tol = inf;

    bigsint = cell(1,length(freqp));
    bigcost = cell(1,length(freqp));

    options = optimoptions(@fminimax, ...
                           'ConstraintTolerance', 1e-4, ...
                           'StepTolerance', 1e-16, ...
                           'FunctionTolerance', 1e-4, ...
                           'Display', 'iter', ...
                           'MaxFunctionEvaluations', 1e4, ...
                           'MaxIterations', 1e4);

    %% step 2 solve the problem
    while(tol > delta)

        % bigsint
        for i = 1:length(freqp)
            for j = 1:length(coeff)
                temp(j) = sin(j*freqp(i)-0.5*length(coeff)*freqp(i)-0.5*phred(i));
            end
            bigsint{i} = temp;
            clear temp;
        end

        % bigcost
        for i = 1:length(freqp)
            for j = 1:length(coeff)
                temp(j) = cos(j*freqp(i)-0.5*length(coeff)*freqp(i)-0.5*phred(i));
            end
            bigcost{i} = temp;
            clear temp;
        end

        efhandle = @(x)errfun(x,coeff,freqp,phred,bigsint,bigcost);
        nchandle = @(x)stbcon(x);

        [coeff,newerr] = fminimax(efhandle,coeffinit,[],[],[],[],[],[],nchandle,options);

        tol = ((max(abs(2*atan(olderr))) - max(abs(2*atan(newerr)))) ...
               / max(abs(2*atan(olderr))));

        olderr = newerr;

    end

    coeff = [1,coeff];
    disp(coeff);

    %% APPENDiX the experation of error function for minimax algorithm('fminimax')

    % error function WITH UNKNOWN 'x' (anonymous)
    function err = errfun(x,coeff,freqp,phred,bigsint,bigcost)
        for i = 1:length(freqp)
            err(i) = (-sin(0.5*(length(x)*freqp(i)+phred(i)))+bigsint{i}*[x]') ...
                     / (abs(cos(0.5*(length(coeff)*freqp(i)+phred(i)))+bigcost{i}*[coeff]'));
        end
    end


    % non-linear constraint for test the stability
    function [c,ceq] = stbcon(x)
        c = [];
        ceq = isstable(fliplr([1,x]),[1,x]) - 1;
    end

end
