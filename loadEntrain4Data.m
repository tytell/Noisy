function data = loadEntrain4Data(filename, varargin)

opt.checksweep = true;
opt.stimuluslocation = [];
opt = parsevarargin(opt,varargin, 4);

if ~exist('TDMS_getStruct','file')
    error('Cannot find TDMS library');
end

fileinfo = TDMS_getStruct(filename,[],{'GET_DATA_OPTION','getnone'});

startdatevec = TDMS_readChannelOrGroup(filename,'Untitled','StartDateVector');
sampfreq = fileinfo.Input.Props.SampFreqHz;
sampfreq = double(sampfreq);
motorfreq = fileinfo.Output.Props.MotorOutputFreqHz;

if ~isfield(fileinfo.Output.Motor.Props,'StimulusType')
    error('Stimulus type info not in file... Weird.');
end

stimtype = fileinfo.Output.Motor.Props.StimulusType;

[sig,channelnames] = TDMS_readChannelOrGroup(filename,'Input');
sig = cat(1,sig{:});
sig = sig';

ang = TDMS_readChannelOrGroup(filename,'Output','Motor');
ang = ang';
if ~isempty(opt.stimuluslocation)
    stimpos = opt.stimuluslocation;
elseif isfield(fileinfo.Output.Motor.Props,'StimulusLocation')
    stimpos = fileinfo.Output.Motor.Props.StimulusLocation;
else
    stimpos = input('Stimulus location? ');
end

t = (0:size(sig,1)-1)' / sampfreq;
dt = 1/sampfreq;

if motorfreq > 100
    motorsavefreq = diground(length(ang)/t(end),10);
    fprintf('Assuming motor output data was saved at %dHz\n', motorsavefreq);
else
    motorsavefreq = motorfreq;
end

tang = (0:size(ang,1)-1)' / motorsavefreq;
ang = interp1(tang,ang, t);

data.sig = sig;
data.t = t;
data.ang = ang;
data.channelnames = channelnames;
data.stimuluspos = stimpos;

if isfield(fileinfo.Output,'SinePhase')
    phase = TDMS_readChannelOrGroup(filename,'Output','SinePhase');
    phase = phase';
    phase = unmod(phase,1);
    
    phase = interp1(tang,phase, t);
    
    cycle = floor(phase);
    phase = mod(phase,1);
elseif isfield(fileinfo.Output,'RampPhase')
    rampphase = TDMS_readChannelOrGroup(filename,'Output','RampPhase');
    rampphase = rampphase';
    
    rampphase = interp1(tang,rampphase, t, 'nearest');
else
    phase = [];
    cycle = [];
end

if isfield(fileinfo.Output,'StimPhase')
    stimphase = TDMS_readChannelOrGroup(filename,'Output','StimPhase');
    stimphase = stimphase';
    
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

stimattrs = fileinfo.Output.Motor.Props;
[stimattrdata,stimattrnames] = TDMS_readChannelOrGroup(filename,'Output');
isattr = ~ismember(stimattrnames,{'Motor','SinePhase','StimPhase','RampPhase'});
if any(isattr)
    attrlen = cellfun(@length,stimattrdata(isattr));
    if (any(attrlen ~= attrlen(1)))
        error('Something is weird with the attributes');
    end
    stimattrdata = stimattrdata(isattr);
    for i = 1:length(stimattrdata)
        stimattrdata{i} = stimattrdata{i}';
    end
    stimattr2 = cell2struct(stimattrdata',stimattrnames(isattr));
    stimattrs = joinstructfields(stimattrs,stimattr2);
end

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
        if isempty(phase)
            stimper = 1/data.SinFreqStartHz;
            
            phase = NaN(size(data.t));
            issig = (data.t >= data.BeforeDurSec) & (data.t < data.StimDurSec-data.AfterDurSec);
            phase(issig) = (data.t(issig) - data.BeforeDurSec)/stimper;
            
            dt = data.t(2) - data.t(1);
            
            indshift = round((data.Time + data.BeforeDurSec + data.Phase*stimper)/dt) + 1;
            indzero = round((data.Time + data.BeforeDurSec + data.Phase*stimper + (1 - data.Phase-data.PhaseChange)*stimper)/dt) + 1;
            
            dphase = zeros(size(phase));
            dphase(indshift) = data.PhaseChange;
            dphase = cumsum(dphase);
            phase = phase + dphase;

            stimphase = NaN(size(data.t));
            stimcycle = NaN(size(data.t));
            tt = [data.Time; data.Time(end)+10*stimper] + data.BeforeDurSec;
            for i = 1:length(tt)-1
                isstim = (data.t >= tt(i)) & (data.t < tt(i+1));
                stimphase(isstim) = phase(isstim);
                stimphase(isstim) = stimphase(isstim) - phase(indzero(i));
                stimcycle(isstim) = i;
            end

            cycle = floor(phase);
            phase = mod(phase,1);
            
                
        end
            
        data.phase = phase;
        data.cycle = cycle;
        data.stimphase = stimphase;
        data.stimcycle = stimcycle;
        
        data.stimfreq = data.SinFreqStartHz;
        data.amp = data.SinAmpStartDeg;
        data.noise = data.NoiseAmpStartDeg;
        
    case 'Ramps'
        data.rampphase = rampphase;
        
    otherwise
        error('Unrecognized treatment type: %s',stimtype);
end
        
        
    

