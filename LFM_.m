%% Ambiguity waveform:

lfmwaveform = phased.LinearFMWaveform('SampleRate',16,'SweepBandwidth',4,'PRF',1,'PulseWidth', 0.25);

release(lfmwaveform);

lfmwaveform.NumPulses = 4;

wav = lfmwaveform();

[afmag,delay,doppler] = ambgfun(wav,lfmwaveform.SampleRate,lfmwaveform.PRF);

%ambgfun(wav,lfmwaveform.SampleRate,lfmwaveform.PRF);


surf(delay, doppler, (afmag));

xlabel('Delay', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Doppler', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('\chi', 'FontSize', 12, 'FontWeight', 'bold');

