%% Equiripple Group Delay Error Allpass Filters Designing
%  log:
%   version1:
%    july 10, 2017
%
%  author:
%   matt ma @scie, shanghai university
%   mattma9209@aliyun.com
%
%  input vars:
%   freqp: frequency points;
%   maxpr: maximum pole radius(0.99 by default);
%   phred: phase response desired.
%
%  output vars:
%   coeff: coefficient vector of designed allpass filter as desired.
%
%  reference:
%   Design and Application of Allpass Filters with Equiripple Group Delay Errors(2013);

function coeff = eqrpgdr(freqp, maxpr, phred)

%% step 1 initialize
coeff = 0;
k = 0;
olderr = inf;

coeffinit = 0; % correct?
lowerbd = -2;
upperbd = +2;
efhandle = @(x)errfun(x,coeff,freqp,phred,bigsint,bigcost);

%% step 2 solve the problem
bigsint = cell(1,length(freqp));
bigcost = cell(1,length(freqp));
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

    [newcoeff,newerr] = fminimax(efhandle,coeffinit,[],[],[],[],-2,+2,);

    if(((max(newerr)-max(olderr))/max(olderr))<delta)
        break;
    end

    olderr = newerr;
    coeff = [coeff,newcoeff];
end

%% APPENDiX the experation of error function for minimax algorithm('fminimax')

% error function WITH UNKNOWN 'x' (anonymous)
    for i = 1:length(freqp)
        err(i) = (-sin(0.5*(length(coeff)*freqp(i)+phred(i)))+bigsint{i}*[coeff,0]') ...
                 / (abs(cos(0.5*(length(coeff)*freqp(i)+phred(i)))+bigcost{i}*[coeff,0]'));
    end
end
