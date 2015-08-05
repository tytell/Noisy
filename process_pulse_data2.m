function data = process_pulse_data2(data, varargin)

opt.sampfreqlo = 50;
opt.nrand = 2;
opt.smoothdur = 0.2;        %sec
opt = parsevarargin(opt,varargin,2);

%get the pulse center times
pulset = data.Time + data.BeforeDurSec + data.Phase/data.stimfreq;
pulsedur = data.Duration(1);

nchan = size(data.sig,2);
sampfreq = 1./(data.t(2) - data.t(1));

%smooth the raw burst recordings
%and downsample
tlo = (0:1/opt.sampfreqlo:data.t(end))';
siglo = NaN(length(tlo),nchan);
for c = 1:nchan
    sigsmooth1 = smooth(data.t,abs(data.sig(:,c))-nanmean(abs(data.sig(:,c))),...
        opt.smoothdur*sampfreq,'loess');
    siglo(:,c) = interp1(data.t,sigsmooth1, tlo);
end

%use phaser to estimate system phase
[~,phi] = newPhaser(siglo');
phi = phi/(2*pi);
freq = deriv(tlo,phi);

%get burst phases
burstphase = interp1(tlo,phi, data.burstt);
[burstnum,burstpulse] = get_burst_pulse_num(data.burston,data.burstoff, ...
    pulset, pulsedur);

%get the pulse on and off phase
pulseonph = interp1(tlo,phi, pulset-pulsedur/2);
pulseoffph = interp1(tlo,phi, pulset+pulsedur/2);

%change in phase during the pulse
dphase = pulseoffph - pulseonph;

%get the frequency after pulses
%one cycle after the pulse
ph1 = pulseoffph + 1;
t1 = interp1(phi,tlo, ph1);
freq1 = 1./(t1 - (pulset+pulsedur/2));

%two cycles after the pulse
ph2 = pulseoffph + 2;
t2 = interp1(phi,tlo, ph2);
freq2 = 1./(t2 - t1);

%get the phase at the center of the pulse
%extrapolate the frequency from the 3 sec before the pulse
pulsephase = NaN(size(pulset));
prefreq = NaN(size(pulset));
for i = 1:length(pulset)
    ispre = (tlo > pulset(i)-pulsedur/2 - 3) & (tlo < pulset(i)-pulsedur/2);
    prefreq(i) = nanmean(freq(ispre));
    
    pulsephase(i) = pulseonph(i) + prefreq(i)*pulsedur/2;
end

%now get a bunch of random times that don't overlap with pulses
%and see how much phase tends to change over that time
randt = [];
while length(randt) < opt.nrand*length(pulset)
    newrand = opt.nrand*length(pulset) - length(randt);
    
    randt1 = rand(newrand,1) * (pulset(end) - pulset(1)) + pulset(1);
    good = true(size(randt1));
    for i = 1:length(randt1)
        m = min(abs(randt1(i) - pulset));
        good(i) = m > pulsedur;
    end
    
    randt = cat(1,randt,randt1(good));
end

randonph = interp1(tlo,phi, randt-pulsedur/2);
randoffph = interp1(tlo,phi, randt+pulsedur/2);
randphase = interp1(tlo,phi, randt);

data = rmfield(data,{'stimcyclet','burstspercycle','burstcycle','burstcyclet'});

data.burstphase = burstphase;
data.burstnum = burstnum;
data.burstpulse = burstpulse;

data.phase = phi;
data.freq = freq;

data.pulseonph = pulseonph;
data.pulseoffph = pulseoffph;
data.pulseph = pulsephase;
data.freqduring = dphase/pulsedur;
data.freq0 = prefreq;
data.freq1 = freq1;
data.freq2 = freq2;

data.randphase = randphase;
data.freqduringrand = (randoffph - randonph) / pulsedur;





