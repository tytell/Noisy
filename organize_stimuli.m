function data = organize_stimuli(data, varargin)

opt.prestimcycles = 5;
opt.basephaselag = 0.01;
opt = parsevarargin(opt, varargin, 2);

nchan = size(data.sig,2);
stimper = 1/data.SinFreqStartHz;
goodchan = data.goodchan > 0;
goodind = find(goodchan);

side = regexp(data.channelnames,'[LR]','once','match');
isleft = cellfun(@(x) (~isempty(x)) && (x == 'L'), side);

seg = regexp(data.channelnames,'\d+','once','match');
seg = cellfun(@str2double,seg);

relseg = data.stimuluspos;

burstrelphase = (seg - relseg)*opt.basephaselag;
burstrelphase(~isleft) = burstrelphase(~isleft) + 0.5;
burstrelphase = burstrelphase(:)';
data.burstrelphase = burstrelphase;

spikestimind = NaN(size(data.spiket));
burststimind = NaN(size(data.burstt));
for c = 1:nchan
    good = isfinite(data.spiket(:,c));
    spikestimind(good,c) = interp1(data.t,data.stimcycle, data.spiket(good,c));
    
    good = isfinite(data.burstt(:,c));
    burststimind(good,c) = interp1(data.t,data.stimcycle, data.burstt(good,c));
end

burston = catuneven(1,data.burst.on,'keepempty')';
burstoff = catuneven(1,data.burst.off,'keepempty')';

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

        t0 = tstart + data.Phase/data.SinFreqStartHz;
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

dt = data.t(2) - data.t(1);
t2 = (0:dt:stimdur)' - opt.prestimcycles/data.SinFreqStartHz;
t2 = t2 - first(t2,t2 >= 0);

data.tstim = t2;
data.sigstim = [];
data.angstim = [];
data.spiketstim = [];
data.bursttstim = [];
data.burstcyclestim = [];
data.burstonstim = [];
data.burstdurstim = [];
data.burstphasestim = [];
data.burstindstim = [];
data.burstprepoststim = [];
data.burstperstim = [];
data.burstphaseatstim = NaN(nstim,1);

for i = 1:nstim
    isstim = (data.t > tstart(i)) & (data.t <= tend(i));
    
    sig1 = data.sig(isstim,:);
    ang1 = data.ang(isstim);
    t1 = data.t(isstim) - t0(i);

    ang2 = interp1(t1,ang1, t2);
    
    spiket1 = [];
    burstt1 = [];
    burston1 = [];
    burstdur1 = [];
    burstphase1 = [];
    burstind1 = [];
    prepost1 = [];
    sig2 = zeros(size(t2,1),nchan);
    for j = 1:nchan
        sig2(:,j) = interp1(t1,sig1(:,j), t2);
        
        isstim1 = (data.spiket(:,j) >= tstart(i)) & (data.spiket(:,j) <= tend(i));
        if any(isstim1)
            spiket1 = catuneven(2,spiket1,data.spiket(isstim1,j));
        else
            spiket1 = catuneven(2,spiket1,NaN);
        end
        
        isstim1 = (data.burstt(:,j) >= tstart(i)) & (data.burstt(:,j) <= tend(i));
        if any(isstim1)
            bon1 = burston(isstim1,j);
            boff1 = burstoff(isstim1,j);

            prepost2 = zeros(size(bon1));
            prepost2(boff1 < t0(i)-data.Duration(i)/2) = -1; 
            prepost2(bon1 > t0(i)+data.Duration(i)/2) = 1;
            k = first(prepost2 == 1);
            if (k+1 <= length(prepost2))
                prepost2(k+1) = 2;
            end
            if (k+2 <= length(prepost2))
                prepost2(k+2:end) = 3;
            end
            
            prepost1 = catuneven(2,prepost1,prepost2);
            
            burstt1 = catuneven(2,burstt1,data.burstt(isstim1,j));
            burston1 = catuneven(2,burston1,bon1);
            burstdur1 = catuneven(2,burstdur1,boff1 - bon1);
            if (isfield(data,'burstphase'))
                burstphase1 = catuneven(2,burstphase1,data.burstphase(isstim1,j));
            else
                burstphase1 = catuneven(2,burstphase1,NaN(sum(isstim1),1));
            end
            burstind1 = catuneven(2,burstind1,find(isstim1));
        else
            burstt1 = catuneven(2,burstt1,NaN);
            burston1 = catuneven(2,burston1,NaN);
            burstdur1 = catuneven(2,burstdur1,NaN);
            burstphase1 = catuneven(2,burstphase1,NaN);
            burstind1 = catuneven(2,burstind1,NaN);
            prepost1 = catuneven(2,prepost1,NaN);
        end
    end
    
    spiket1 = spiket1 - t0(i);
    burstt1 = burstt1 - t0(i);
    burston1 = burston1 - t0(i);
    
    burstper1 = NaN(size(burstt1));
    burstper1(2:end,:) = diff(burstt1);
    
    lastper = burstper1;
    lastper(prepost1 ~= -1) = NaN;
    lastper = nanmedian(flatten(lastper(:,goodchan)));
    
    stimphase1 = (0 - last(burston1,prepost1 == -1)) / lastper;
    stimphase1 = stimphase1 + burstrelphase;
    stimphase1 = angmean(2*pi*stimphase1(goodchan)) / (2*pi);
    
    if all(isnan(burstphase1(:)))
        burstphase1 = burston1 / lastper + stimphase1;
    else
        burstphase1 = burston1 / stimper + stimphase1;
    end
    
    if ~issamedur
        data.tstim = catuneven(3,data.tstim,t2);
    end
    data.sigstim = catuneven(3,data.sigstim,sig2);
    data.angstim = catuneven(3,data.angstim,ang2);
    data.spiketstim = catuneven(3,data.spiketstim,spiket1);
    data.bursttstim = catuneven(3,data.bursttstim,burstt1);
    data.burstcyclestim = catuneven(3,data.burstcyclestim, floor(burstt1/stimper));
    data.burstonstim = catuneven(3,data.burstonstim,burston1);
    data.burstdurstim = catuneven(3,data.burstdurstim,burstdur1);
    data.burstphasestim = catuneven(3,data.burstphasestim,burstphase1);
    data.burstprepoststim = catuneven(3,data.burstprepoststim,prepost1);
    data.burstindstim = catuneven(3,data.burstindstim,burstind1);
    data.burstperstim = catuneven(3,data.burstperstim,burstper1);
    data.burstphaseatstim(i) = stimphase1;
end


    
    
    