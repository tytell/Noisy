function data = loadEntrainData(filename,treatments,ind, varargin)

opt.output = 'motor';       % or encoder

opt = parsevarargin(opt,varargin, 4);

[startdatevec,err] = hdf5err(@h5readatt,filename,'/','StartDateVector');
[sampfreq,err] = hdf5err(@h5readatt,filename,'/Input','SampFreqHz');
sampfreq = double(sampfreq);
[updatefreq,err] = hdf5err(@h5readatt,filename,'/Input','UpdateFreqHz');
updatefreq = double(updatefreq);

chunksz = sampfreq / updatefreq;

startsamp = treatments.startsample(ind);
if (ind < length(treatments.startsample))
    nsamp = treatments.startsample(ind+1) - startsamp - 1;
else
    nsamp = Inf;
end

sig = h5read(filename,'/Input/Voltage', [startsamp 1], ...
    [nsamp Inf]);

switch opt.output
    case 'motor'
        ang = h5read(filename,'/Output/Motor', [startsamp+3*chunksz 2], ...
            [nsamp 1]);
        
    case 'encoder'
        ang = h5read(filename,'/Input/Encoder', treatments.startsample(ind), ...
            nsamp);
end

t = startsamp/sampfreq + (0:size(sig,1)-1)' / sampfreq;

data.sig = sig;
data.t = t;
data.ang = ang;

if (strcmp(treatments.type{ind}, 'Sine'))
    phase = mod((t - startsamp) * treatments.frequencyhz(ind), 1);
    data.phase = phase;
    data.stimfreq = treatments.frequencyhz(ind);
end


        
        
    

