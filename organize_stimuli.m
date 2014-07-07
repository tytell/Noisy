function data = organize_stimuli(data, varargin)

opt.prestimcycles = 5;
opt = parsevarargin(opt, varargin, 2);

nchan = size(data.sig,2);
stimper = 1/data.SinFreqStartHz;

spikestimind = NaN(size(data.spiket));
burststimind = NaN(size(data.burstt));
for i = 1:nchan
    good = isfinite(data.spiket(:,i));
    spikestimind(good,i) = interp1(data.t,data.stimcycle, data.spiket(good,i));
    
    good = isfinite(data.burstt(:,i));
    burststimind(good,i) = interp1(data.t,data.stimcycle, data.burstt(good,i));
end

burston = catuneven(1,data.burst.on)';
burstoff = catuneven(1,data.burst.off)';

stimdur = mode(round(diff(data.Time)/stimper)) * stimper;

nstim = length(data.Time);
switch data.StimulusType
    case 'Pulses'
        %check that we're close to an even number of cycles, then round
        %to an even number of cycles
        tstart = data.Time;
        assert(all(abs(round(tstart/stimper)*stimper - tstart) < 0.1*stimper));
        tstart = round(tstart/stimper) * stimper + data.BeforeDurSec;
        tend = tstart + stimdur;

        t0 = tstart + 1/data.SinFreqStartHz;    % 1 period later
        tstart = tstart - opt.prestimcycles/data.SinFreqStartHz;
        issamedur = true;

    case 'Shifts'
        stimdur = stimdur - stimper;
        
        t0 = data.Time + data.BeforeDurSec + data.Phase*stimper + ...
            (1 - data.Phase-data.PhaseChange)*stimper;
        tend = t0 + stimdur;
        tstart = t0 - data.PhaseChange*stimper - opt.prestimcycles*stimper;
        issamedur = false;
end

data.tstim = [];
data.sigstim = [];
data.angstim = [];
data.spiketstim = [];
data.bursttstim = [];
data.burstcyclestim = [];
data.burstonstim = [];
data.burstdurstim = [];
data.burstphasestim = [];
data.burstindstim = [];
for i = 1:nstim
    isstim = (data.t > tstart(i)) & (data.t <= tend(i));
    
    sig1 = data.sig(isstim,:);
    ang1 = data.ang(isstim);
    t1 = data.t(isstim) - t0(i);
    
    spiket1 = [];
    burstt1 = [];
    burston1 = [];
    burstdur1 = [];
    burstphase1 = [];
    burstind1 = [];
    for j = 1:nchan
        isstim1 = (data.spiket(:,j) >= tstart(i)) & (data.spiket(:,j) <= tend(i));
        spiket1 = catuneven(2,spiket1,data.spiket(isstim1,j));
        
        isstim1 = (data.burstt(:,j) >= tstart(i)) & (data.burstt(:,j) <= tend(i));
        if any(isstim1)
            burstt1 = catuneven(2,burstt1,data.burstt(isstim1,j));
            burston1 = catuneven(2,burston1,burston(isstim1,j));
            burstdur1 = catuneven(2,burstdur1,burstoff(isstim1,j) - burston(isstim1,j));
            burstphase1 = catuneven(2,burstphase1,data.burstphase(isstim1,j));
            burstind1 = catuneven(2,burstind1,find(isstim1));
        else
            burstt1 = catuneven(2,burstt1,NaN);
            burston1 = catuneven(2,burston1,NaN);
            burstdur1 = catuneven(2,burstdur1,NaN);
            burstphase1 = catuneven(2,burstphase1,NaN);
            burstind1 = catuneven(2,burstind1,NaN);
        end
    end
    
    spiket1 = spiket1 - t0(i);
    burstt1 = burstt1 - t0(i);
    burston1 = burston1 - t0(i);
    
    if issamedur
        data.tstim = t1;
    else
        data.tstim = catuneven(3,data.tstim,t1);
    end
    data.sigstim = catuneven(3,data.sigstim,sig1);
    data.angstim = catuneven(3,data.angstim,ang1);
    data.spiketstim = catuneven(3,data.spiketstim,spiket1);
    data.bursttstim = catuneven(3,data.bursttstim,burstt1);
    data.burstcyclestim = catuneven(3,data.burstcyclestim, floor(burstt1));
    data.burstonstim = catuneven(3,data.burstonstim,burston1);
    data.burstdurstim = catuneven(3,data.burstdurstim,burstdur1);
    data.burstphasestim = catuneven(3,data.burstphasestim,burstphase1);
    data.burstindstim = catuneven(3,data.burstindstim,burstind1);
end


    
    
    