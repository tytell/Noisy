function data = loadEntrain2Data(filename, varargin)

opt.checksweep = true;

opt = parsevarargin(opt,varargin, 4);

[startdatevec,err] = hdf5err(@h5readatt,filename,'/','StartDateVector');
[sampfreq,err] = hdf5err(@h5readatt,filename,'/Input','SampFreqHz');
sampfreq = double(sampfreq);
[motorfreq,err] = hdf5err(@h5readatt,filename,'/Output','MotorOutputFreqHz');
motorfreq = double(motorfreq);

sig = h5read(filename,'/Input/Voltage');
ang = h5read(filename,'/Output/Motor');

t = (0:size(sig,1)-1)' / sampfreq;
dt = 1/sampfreq;

tang = (0:size(ang,1)-1)' / motorfreq;
ang = interp1(tang,ang, t);

data.sig = sig;
data.t = t;
data.ang = ang;

[stimtype,err] = hdf5err(@h5readatt,filename,'/Output','StimulusType');
stimtype = stimtype{1};

stiminfo = h5info(filename, '/Output');

if ismember('SinePhase',{stiminfo.Datasets.Name})
    phase = h5read(filename,'/Output/SinePhase');
    phase = unmod(phase,1);
    
    phase = interp1(tang,phase, t);
    
    cycle = floor(phase);
    phase = mod(phase,1);
else
    phase = [];
    cycle = [];
end

if ismember('StimPhase',{stiminfo.Datasets.Name})
    stimphase = h5read(filename,'/Output/StimPhase');
    iswrap = diff(stimphase) < 0;
    wrapval = mode(ceil(stimphase(iswrap)));
    
    stimphase = unmod(stimphase,wrapval);
    stimphase = interp1(tang,stimphase, t);
    
    stimcycle = floor(stimphase/wrapval);
    stimphase = mod(stimphase,wrapval);
else
    stimphase = [];
    stimcycle = [];
end

stimattrs = cell2struct({stiminfo.Attributes.Value}',{stiminfo.Attributes.Name}');
data = joinstructfields(data,stimattrs);
data.StimulusType = stimtype;

t = t - data.BeforeDurSec;
isstim = (t >= 0) & (t < data.StimDurSec);
n = sum(isstim);

switch stimtype
    case 'Zero'
        data.phase = NaN(size(t));
        data.cycle = zeros(size(t));
        data.stimfreq = NaN;
        data.amp = 0;
        data.noise = 0;
        
    case 'Sine'
        if ~isempty(phase)
            data.phase = phase;
            data.cycle = cycle;
        else
            data.phase = t*data.SinFreqStartHz;
            data.cycle = floor(data.phase)+1;
            data.phase = mod(data.phase,1);
        end
        data.stimfreq = data.SinFreqStartHz;
        data.amp = data.SinAmpStartDeg;
        data.noise = data.NoiseAmpStartDeg;
        data.stimphase = data.phase;
        data.stimcycle = ones(size(data.phase));
        
    case 'FrequencySweep'
        if isempty(phase)
            %Chirp Pattern Details
            %If the sequence Y represents Chirp Pattern, the Chirp Pattern VI
            %obtains the elements of Y using the following equation:
            %
            %yi = A*sin((0.5*a*i + b)*i)
            %
            %for i = 0, 1, 2, …, n – 1
            %
            % where A is amplitude,
            % 
            % a = 2(f2 – f1)/n,
            % 
            % b = 2f1,
            % 
            % f1 is the beginning frequency in normalized units of cycles/sample,
            % 
            % f2 is the ending frequency in normalized units of cycles/sample,
            % 
            % n is the number of samples.

            f1 = data.SinFreqStartHz * dt;
            f2 = data.SinFreqEndHz * dt;
            k = (0:n-1)';

            a = 2*(f2 - f1)/n;
            b = 2*f1;

            phase = NaN(size(t));
            phase(isstim) = 0.5*(0.5*a*k + b).*k;
        
            data.phase = mod(phase,1);
            data.cycle = floor(phase)+1;
        else
            data.phase = phase;
            data.cycle = cycle;
            data.stimphase = stimphase;
            data.stimcycle = stimcycle;
        end
        data.stimfreq = deriv(t,phase);
        data.amp = data.SinAmpStartDeg;
        data.noise = data.NoiseAmpStartDeg;
    
    case 'AmplitudeSweep'
        if isempty(phase)
            phase = NaN(size(t));
            phase(isstim) = t(isstim) * data.SinFreqStartHz;
            data.phase = mod(phase,1);
            data.cycle = floor(phase)+1;
        else
            data.phase = phase;
            data.cycle = cycle;
            data.stimphase = stimphase;
            data.stimcycle = stimcycle;
        end            
        data.stimfreq = data.SinFreqStartHz;
        
        if (data.SinAmpStartDeg == data.SinAmpEndDeg)
            data.amp = data.SinAmpStartDeg;
        else
            data.amp = zero(size(t));
            data.amp(isstim) = t(isstim)*...
                (data.SinAmpEndDeg-data.SinAmpStartDeg)/data.StimDurSec + ...
                data.SinAmpStartDeg;
        end
        
        if (data.NoiseAmpStartDeg == data.NoiseAmpEndDeg)
            data.noise = data.NoiseAmpStartDeg;
        else
            data.noise = zeros(size(t));
            data.noise(isstim) = t(isstim)*...
                (data.NoiseAmpEndDeg-data.NoiseAmpStartDeg)/data.StimDurSec + ...
                data.NoiseAmpStartDeg;
        end
        
    case 'Pulses'
        data.phase = phase;
        data.stimphase = stimphase;
        data.cycle = cycle;
        data.stimcycle = stimcycle;

        if ceil(max(stimphase)) == 1
            cyclesperstim = mode(ceil(diff(data.Time)*data.SinFreqStartHz));
            data.stimphase = data.stimphase * cyclesperstim - 1;
        end
        
        data.stimfreq = data.SinFreqStartHz;
        data.amp = data.SinAmpStartDeg;
        data.noise = data.NoiseAmpStartDeg;

    case 'Shifts'
        data.phase = phase;
        data.stimphase = stimphase;
        data.cycle = cycle;
        data.stimcycle = stimcycle;

        if ceil(max(stimphase)) == 1
            cyclesperstim = mode(ceil(diff(data.Time)*data.SinFreqStartHz));
            data.stimphase = data.stimphase * cyclesperstim - 1;
        end
        
        data.stimfreq = data.SinFreqStartHz;
        data.amp = data.SinAmpStartDeg;
        data.noise = data.NoiseAmpStartDeg;
    otherwise
        error('Unrecognized treatment type: %s',stimtype);
end
        
        
    

