filename = 'L16-caud023.h5';
treatments = loadEntrainTreatments(filename);

data = loadEntrainData(filename,treatments,1);
data = findbursts_gui(data, 'threshold', [0.410379 0.225777 0.138209 ], ...
    'interburstdur', [0.200000 0.200000 0.200000 ], 'minspikes', [3.000000 3.000000 2.000000 ], 'quiet');

[mn1,R1] = angmean(2*pi*data.burstphase);
phasemn(1,:) = mod(mn1 / (2*pi), 1);
R(1,:) = R1;

