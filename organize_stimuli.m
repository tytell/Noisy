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

nstim = max(data.stimcycle);
data.tstim = [];
data.sigstim = [];
data.angstim = [];
data.spiketstim = [];
data.bursttstim = [];
data.burstcyclestim = [];
data.burstonstim = [];
data.burstdurstim = [];
data.burstphasestim = [];
for i = 1:nstim
    tstart = data.Time(i);
    %check that we're close to an even number of cycles, then round
    %to an even number of cycles
    assert(abs(round(tstart/stimper)*stimper - tstart) < 0.1*stimper);
    tstart = round(tstart/stimper) * stimper + data.BeforeDurSec;
    tend = tstart + stimdur;
    
    t0 = tstart + 1/data.SinFreqStartHz;    % 1 period later
    tstart = tstart - opt.prestimcycles/data.SinFreqStartHz;
    
    isstim = (data.t > tstart) & (data.t <= tend);
    
    sig1 = data.sig(isstim,:);
    ang1 = data.ang(isstim);
    t1 = data.t(isstim) - t0;
    
    spiket1 = [];
    burstt1 = [];
    burston1 = [];
    burstdur1 = [];
    burstphase1 = [];
    for j = 1:nchan
        isstim1 = (data.spiket(:,j) >= tstart) & (data.spiket(:,j) <= tend);
        spiket1 = catuneven(2,spiket1,data.spiket(isstim1,j));
        
        isstim1 = (data.burstt(:,j) >= tstart) & (data.burstt(:,j) <= tend);
        if any(isstim1)
            burstt1 = catuneven(2,burstt1,data.burstt(isstim1,j));
            burston1 = catuneven(2,burston1,burston(isstim1,j));
            burstdur1 = catuneven(2,burstdur1,burstoff(isstim1,j) - burston(isstim1,j));
            burstphase1 = catuneven(2,burstphase1,data.burstphase(isstim1,j));
        else
            burstt1 = catuneven(2,burstt1,NaN);
            burston1 = catuneven(2,burston1,NaN);
            burstdur1 = catuneven(2,burstdur1,NaN);
            burstphase1 = catuneven(2,burstphase1,NaN);
        end
    end
    
    spiket1 = spiket1 - t0;
    burstt1 = burstt1 - t0;
    burston1 = burston1 - t0;
    
    data.tstim = t1;
    data.sigstim = catuneven(3,data.sigstim,sig1);
    data.angstim = catuneven(3,data.angstim,ang1);
    data.spiketstim = catuneven(3,data.spiketstim,spiket1);
    data.bursttstim = catuneven(3,data.bursttstim,burstt1);
    data.burstcyclestim = catuneven(3,data.burstcyclestim, floor(burstt1));
    data.burstonstim = catuneven(3,data.burstonstim,burston1);
    data.burstdurstim = catuneven(3,data.burstdurstim,burstdur1);
    data.burstphasestim = catuneven(3,data.burstphasestim,burstphase1);
end


    
    
    