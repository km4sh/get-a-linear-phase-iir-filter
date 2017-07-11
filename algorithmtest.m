% data for test initialize
coeff = rand(1,6);
freqp = [0,0.1,0.3,0.5,0.7,0.9,1];
phred = [-3,-5,-7,-5,-7,-9,-3];

% procedure fot test 
N = length(coeff);
w = freqp;
thetad = phred;
betad = 0.5*(length(coeff).*freqp+phred);

bigsint = cell(1,length(freqp));
bigcost = cell(1,length(freqp));

% bigsint
for m = 1:length(freqp)
    for n = 1:length(coeff)+1
        temp(n) = sin(n*freqp(m)-0.5*length(coeff)*freqp(m)-0.5*phred(m));
    end
    bigsint{m} = temp;
    clear temp;
end

% bigcost
for m = 1:length(freqp)
    for n = 1:length(coeff)+1
        temp(n) = cos(n*freqp(m)-0.5*length(coeff)*freqp(m)-0.5*phred(m));
    end
    bigcost{m} = temp;
    clear temp;
end

% error function
for m = 1:length(freqp)
    err(m) = (-sin(0.5*(length(coeff)*freqp(m)+phred(m)))+bigsint{m}*[coeff,0]') ...
             / (abs(cos(0.5*(length(coeff)*freqp(m)+phred(m)))+bigcost{m}*[coeff,0]'));
end

% display
err


%% just take a room to write the constraint for isstable
%  #### nonlcon should like this
function [c,ceq] = stbcon(x,coeff)
    c = [];
    ceq = isstable(fliplr([coeff,x]),[coeff,x]) - 1;
end
