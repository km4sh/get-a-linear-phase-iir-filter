%% example 1

freqp = (0:320)*pi/400;

phred = -7.0516*freqp;

coeff = eqrpgdr(freqp,phred, 8);
