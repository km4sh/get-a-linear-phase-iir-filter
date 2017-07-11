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

%cutrange = min(min(winrange1),min(winrange2)):max(max(winrange1),max(winrange2));



for i = 1:32
    for j = 0:2048:65536-512
        cutrange = 1:max(max(winrange1),max(winrange2))+j;

        tempdevbrir1 = devbrir1(cutrange);
        tempdevbrir2 = devbrir2(cutrange);

        freqr1 = fft(devbrir1);
        freqr2 = fft(devbrir2);
        phasr1 = phase(freqr1);
        phasr2 = phase(freqr2);
    end
    phasr(:,:,i) = [phasr1',phasr2'];
end
