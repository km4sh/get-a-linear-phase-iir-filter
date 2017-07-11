load rir.mat;

[max1,maxp1] = max(rir1.^2);
[max2,maxp2] = max(rir2.^2);

hannwin = hann(256);

devbrir1 = ones(1,length(rir1))*hannwin(1);
devbrir2 = ones(1,length(rir2))*hannwin(1);

winrange1 = maxp1-ceil(length(hannwin)/2):maxp1+ ...
    floor(length(hannwin)/2)-1;
winrange2 = maxp2-ceil(length(hannwin)/2):maxp2+ ...
    floor(length(hannwin)/2)-1;

devbrir1(winrange1) = rir1(winrange1).*hannwin;
devbrir2(winrange2) = rir2(winrange2).*hannwin;

i = 1;
phasr = cell(1,25);
for j = 16384:512:28672
    cutrange = 1:max(max(winrange1),max(winrange2))+j;

    tempdevbrir1 = devbrir1(cutrange);
    tempdevbrir2 = devbrir2(cutrange);

    freqr1 = fft(tempdevbrir1);
    freqr2 = fft(tempdevbrir2);
    phasr1 = phase(freqr1);
    phasr2 = phase(freqr2);
    phasr{i} = [phasr1',phasr2'];
    i = i+1;
    i
end
