%% example 1

freqp = (0:320)*pi/400;

phred = -7.0615*freqp;

coeff = eqrpgdr(freqp,phred, 8);