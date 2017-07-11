% data for test initialize
%coeff = rand(1,6);
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

% error function
for i = 1:length(freqp)
    err(i) = (-sin(0.5*(length(coeff)*freqp(i)+phred(i)))+bigsint{i}*[coeff,0]') ...
             / (abs(cos(0.5*(length(coeff)*freqp(i)+phred(i)))+bigcost{i}*[coeff,0]'));
end

% display
err

