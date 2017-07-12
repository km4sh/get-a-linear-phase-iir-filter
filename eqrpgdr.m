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

function [coeff,errbuf] = eqrpgdr(freqp, phred)

    %% step 1 initialize
    coeff = 0;
    olderr = inf;
    delta = 1e-6;
    maxerr = inf;
    errbuf = [];

    bigsint = cell(1,length(freqp));
    bigcost = cell(1,length(freqp));

    options = optimoptions(@fminimax, ...
                           'ConstraintTolerance', 1e-4, ...
                           'StepTolerance', 1e-16, ...
                           'FunctionTolerance', 1e-4, ...
                           'Display', 'off', ...
                           'MaxFunctionEvaluations', 1e4, ...
                           'MaxIterations', 1e4);
tic
    %% step 2 solve the problem
    while(length(coeff)<12)

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

        for x = -1.5:1e-3:1.5

            fprintf('\nprocessing');
            fprintf('%7.3f',[coeff,x]);

            curmaxerr = max(errfun(x,coeff,freqp,phred,bigsint,bigcost));

            fprintf(' #MAXERR curnt = %6.1f std = %6.1f',curmaxerr,maxerr);

            if(length(coeff) == 1)
                stableflag = isstable(fliplr([1,x]),[1,x]);
            else
                stableflag = isstable( ...
                    fliplr([1,coeff(2:length(coeff)),x]),[1,coeff(2:length(coeff)),x]);
            end

            if((curmaxerr < maxerr) && stableflag)
                newcoeff = x;
                maxerr = curmaxerr;
                fprintf(' # refreshed!')
            elseif((curmaxerr > maxerr))
                break;
            end
        end

        newerr = maxerr;
        maxerr = inf;
        errbuf = [errbuf,newerr];
        tol = ((max(abs(2*atan(newerr))) - max(abs(2*atan(olderr)))) ...
               / max(abs(2*atan(olderr))));

        if(tol < delta)
            fprintf('\nstep out by error tolerance never changes!\n');
            fprintf('error tolerance fail with %12.6f\n',tol);
            break;
        end

        olderr = newerr;
        coeff = [coeff,newcoeff];

    end
toc
    a = [1,coeff(2:length(coeff))];
    b = fliplr(a);
    coeff = [b;a];

    %% APPENDiX the experation of error function for minimax algorithm('fminimax')

    % error function WITH UNKNOWN 'x' (anonymous)
    function err = errfun(x,coeff,freqp,phred,bigsint,bigcost)
        for i = 1:length(freqp)
            err(i) = abs(-sin(0.5*(length(coeff)*freqp(i)+phred(i)))+bigsint{i}*[coeff,x]') ...
                     / abs(cos(0.5*(length(coeff)*freqp(i)+phred(i)))+bigcost{i}*[coeff,0]');
        end
    end


    % non-linear constraint for test the stability
    function [c,ceq] = stbcon(x,coeff)
        c = [];
        if(length(coeff) == 1)
            ceq = 0;
        else
            ceq = isstable(fliplr([1,coeff(2:length(coeff)),x]),[1,coeff(2:length(coeff)),x]) - 1;
        end
    end

end


% 1. how to process the first element of 'coeff' (that is 0) well..

% 2. constraint from line 71 to 76 is a problem..

% 3. why fminimax stopped so quickly..

% 4. can the first point of 'freqp' be zero..
