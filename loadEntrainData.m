function data = loadEntrainData(filename,treatments,ind, varargin)

opt.output = 'encoder';       % or encoder
opt.checksweep = true;

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
startt = t(1);
dt = 1/sampfreq;

data.sig = sig;
data.t = t;
data.ang = ang;

switch treatments.type{ind}
    case 'Sine'
        phase = (t - startt) * treatments.frequencyhz(ind);
        data.phase = mod(phase,1);
        data.cycle = floor(phase)+1;
        data.stimfreq = treatments.frequencyhz(ind);
        data.amp = treatments.amplitudedeg(ind);
        data.noise = treatments.noisestddeg(ind);
        
    case 'Frequency Sweep'
        phase = treatments.frequencyhz(ind) * ...
            (treatments.kfreq(ind).^(t-startt) - 1) / log(treatments.kfreq);
        data.phase = mod(phase,1);
        data.cycle = floor(phase)+1;
        
        if (opt.checksweep)
            upind = find((ang(2:end) > 0) & (ang(1:end-1) <= 0));
            downind = find((ang(2:end) <= 0) & (ang(1:end-1) > 0));
            
            downphase = mod(phase(downind),1);
            [~,R] = angmean(2*pi*downphase);
            if (R < 0.9)
                warning('Sweep phase seems to be weird.  Estimating empirically');
                
                %linearly interpolate the true zero crossings
                tup = t(upind) + dt/(ang(upind+1) - ang(upind)) * (0 - ang(upind));
                tdown = t(downind) + dt/(ang(downind+1) - ang(downind)) * (0 - ang(downind));
                
                %ascending zero crossing should have phase 0,1,...
                phup = (0:length(tup)-1)';
                %descending should be phase 0.5,1.5,...
                phdown = (0.5:length(tdown))';
                
                tzero = [tup; tdown];
                [tzero,ord] = sort(tzero);
                phzero = [phup; phdown];
                phzero = phzero(ord);
                phzero = unmod(phzero,1);
                
                span = (t >= tzero(1)) & (t <= tzero(end));
                phase = NaN(size(t));
                phase(span) = interp1(tzero,phzero, t(span));                
            end
        end
        data.phase = mod(phase,1);
        data.cycle = floor(phase)+1;
        data.stimfreq = deriv(t,phase);
        data.amp = treatments.amplitudedeg(ind);
    
    case 'AmplitudeSweep'
        phase = (t - startt) * treatments.frequencyhz(ind);
        data.phase = mod(phase,1);
        data.cycle = floor(phase)+1;
        data.stimfreq = treatments.frequencyhz(ind);
        if (treatments.amplitudechangedegpersec(ind) ~= 0)
            data.amp = treatments.amplitudedeg(ind) + (t - startt) * treatments.amplitudechangedegpersec(ind);
        else
            data.amp = treatments.amplitudedeg(ind);
        end
        if (treatments.noisechangedegpersec(ind) ~= 0)
            data.noise = treatments.noisestddeg(ind) + (t - startt) * treatments.noisechangedegpersec(ind);
        else
            data.noise = treatments.noisestddeg(ind);
        end        
    otherwise
        error('Unrecognized treatment type: %s',treatments.type{ind});
end


        
        
    

