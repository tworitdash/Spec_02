%% (Zero delay response: Doppler cut)

lfmwaveform = phased.LinearFMWaveform('SampleRate',16,'SweepBandwidth',4,'PRF',1,'PulseWidth', 0.25);


release(lfmwaveform);

lfmwaveform.NumPulses = 4;

wav = lfmwaveform();

figure(1);
ambgfun(wav,lfmwaveform.SampleRate,lfmwaveform.PRF, 'Cut', 'Delay');
title('Zero delay response (Doppler Cut)', 'FontSize', 12, 'FontWeight', 'bold')

figure(2);
ambgfun(wav,lfmwaveform.SampleRate,lfmwaveform.PRF, 'Cut', 'Doppler');
title('Zero doppler response (Range Cut)', 'FontSize', 12, 'FontWeight', 'bold');

%xlabel('\tau/T', 'FontSize', 12, 'FontWeight', 'bold');
%ylabel('\chi(0, \nu)', 'FontSize', 12, 'FontWeight', 'bold');
