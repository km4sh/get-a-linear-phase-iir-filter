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

function coeff = eqrpgdr(freqp, phred)

    %% step 1 initialize
    coeff = 0;
    olderr = inf;
    delta = 1e-8;
    coeffinit = zeros(1,8); % correct?

    bigsint = cell(1,length(freqp));
    bigcost = cell(1,length(freqp));

    options = optimoptions(@fminimax, ...
                           'ConstraintTolerance', 1e-4, ...
                           'StepTolerance', 1e-16, ...
                           'FunctionTolerance', 1e-4, ...
                           'Display', 'off', ...
                           'MaxFunctionEvaluations', 1e4, ...
                           'MaxIterations', 1e4);

    %% step 2 solve the problem
    while(length(coeff)<=16)

        % bigsint
        for i = 1:length(freqp)
            for j = 1:length(coeff)+1
                temp(j) = sin(j*freqp(i)-0.5*length(coeff)*freqp(i)-0.5*phred(i));
            end
            bigsint{i} = temp;
            clear temp;
        end

        % bigcost
        for i = 1:length(freqp)
            for j = 1:length(coeff)+1
                temp(j) = cos(j*freqp(i)-0.5*length(coeff)*freqp(i)-0.5*phred(i));
            end
            bigcost{i} = temp;
            clear temp;
        end

        efhandle = @(x)errfun(x,coeff,freqp,phred,bigsint,bigcost);
        nchandle = @(x)stbcon(x,coeff);

        [newcoeff,newerr] = fminimax(efhandle,coeffinit,[],[],[],[],[],[],nchandle,options);

        tol = ((max(abs(2*atan(newerr))) - max(abs(2*atan(olderr)))) ...
               / max(abs(2*atan(olderr))));

        if(tol < delta)
            fprintf('\nstep out by error evaluation never changes!\n');
            fprintf('error tolerance fail with %12.6f\n',tol);
            break;
        end

        olderr = newerr;
        coeff = [coeff,newcoeff];

    end


    %% APPENDiX the experation of error function for minimax algorithm('fminimax')

    % error function WITH UNKNOWN 'x' (anonymous)
    function err = errfun(x,coeff,freqp,phred,bigsint,bigcost)
        for i = 1:length(freqp)
            err(i) = (-sin(0.5*(length(coeff)*freqp(i)+phred(i)))+bigsint{i}*[coeff,x]') ...
                     / (abs(cos(0.5*(length(coeff)*freqp(i)+phred(i)))+bigcost{i}*[coeff,0]'));
        end
    end


    % non-linear constraint for test the stability
    function [c,ceq] = stbcon(x,coeff)
        c = [];
        ceq = isstable(fliplr([coeff,x]),[coeff,x]) - 1;
    end

end
