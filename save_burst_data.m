function save_burst_data(data, outfile)

if (nargin == 1) || isempty(outfile)
    [fn,pn] = uiputfile('*.csv');
    outfile = fullfile(pn,fn);
end

switch data.StimulusType
    case 'Pulses'
        save_pulse_data(data, outfile);
end


function save_pulse_data(data, outfile)

stimper = 1 / data.SinFreqStartHz;

ncycle = mode(flatten(max(data.burstcyclestim)));

onphase = mod(data.burstonstim,stimper);

nchan = size(data.sig,2);
mnphase = zeros(1,nchan);
mnonphase = zeros(1,nchan);
for c = 1:nchan
    burstcycle1 = squeeze(data.burstcyclestim(:,c,:));
    burstphase1 = squeeze(data.burstphasestim(:,c,:));
    goodcyc = (burstcycle1 < -1) | ...
        ((burstcycle1 > 4) & (burstcycle1 <= ncycle-1));
    mnphase1 = angmean(2*pi*burstphase1(goodcyc));
    mnphase(c) = mod(mnphase1/(2*pi),1);
    
    onphase1 = onphase(:,c,:);
    mnonphase1 = angmean(2*pi*onphase1(goodcyc));
    mnonphase(c) = mod(mnonphase1/(2*pi),1);
end

burstfreq = NaN(size(data.bursttstim));
burstfreq(1:end-1,:,:) = 1 ./ diff(data.bursttstim,[],1);

burstint0 = catuneven(1,data.burst.int)';
burstspikerate0 = catuneven(1,data.burst.spikerate)';

burstspikerate = NaN(size(burstfreq));
good = isfinite(data.burstindstim);
burstspikerate(good) = burstspikerate0(data.burstindstim(good));

cycle = data.burstcyclestim;
phase = data.burstphasestim;
prepost = data.burstprepoststim;

expectedt = (cycle + repmat(mnphase,[size(cycle,1) 1 size(cycle,3)]))*stimper;
dphase = (data.bursttstim - expectedt)/stimper;

expectedon = (cycle + repmat(mnonphase,[size(cycle,1) 1 size(cycle,3)]))*stimper;
donphase = (data.burstonstim - expectedon)/stimper;

%look for extra bursts
isextra = cat(1,false(1,nchan,size(cycle,3)), diff(cycle,[],1) == 0);
cycle(isextra) = cycle(isextra) + 0.5;

cycle = permute(cycle,[1 3 2]);
phase = permute(phase,[1 3 2]);
prepost = permute(prepost,[1 3 2]);
burstfreq = permute(burstfreq,[1 3 2]);
dphase = permute(dphase,[1 3 2]);
donphase = permute(donphase,[1 3 2]);
dur = permute(data.burstdurstim / stimper, [1 3 2]);
burstspikerate = permute(burstspikerate, [1 3 2]);

cycle2 = unique(cycle(~isnan(cycle)));
phase2 = NaN(length(cycle2),size(phase,2),nchan);
prepost2 = NaN(size(phase2));
burstfreq2 = NaN(size(phase2));
dphase2 = NaN(size(phase2));
donphase2 = NaN(size(phase2));
dur2 = NaN(size(phase2));
burstspikerate2 = NaN(size(phase2));
for c = 1:nchan
    for i = 1:size(phase,2)
        for j = 1:length(cycle2)
            iscycle = cycle(:,i,c) == cycle2(j);
            if sum(iscycle) == 1
                phase2(j,i,c) = phase(iscycle,i,c);
                prepost2(j,i,c) = prepost(iscycle,i,c);
                burstfreq2(j,i,c) = burstfreq(iscycle,i,c);
                dphase2(j,i,c) = dphase(iscycle,i,c);
                donphase2(j,i,c) = donphase(iscycle,i,c);
                dur2(j,i,c) = dur(iscycle,i,c);
                burstspikerate2(j,i,c) = burstspikerate(iscycle,i,c);
            end
        end
    end
end

cycle2 = repmat(cycle2,[1 size(phase,2)]);

C = cat(4,prepost2,phase2,burstfreq2,dphase2, donphase2, dur2, burstspikerate2);

stimphase = repmat(data.Phase',[size(cycle2,1) 1]);
stimdirec = repmat(data.Direction',[size(cycle2,1) 1]);
stimamp = repmat(data.Amplitude',[size(cycle2,1) 1]);
stimdur = repmat(data.Duration',[size(cycle2,1) 1]);
stimnum = repmat(1:length(data.Phase),[size(cycle2,1) 1]);

X = [stimnum(:) stimphase(:) stimdirec(:) stimamp(:) stimdur(:) cycle2(:)];
lab = {'StimNum','StimPhase','StimDirec','StimAmp','StimDur','CycleNum'};
chanlab = {'CycleType','Phase','Freq','DPhase','DOnPhase','Dur','SpikeRate'};


tplt = '%d,%.3f,%d,%.1f,%.2f,%.1f';
Ctplt = '%d,%.3f,%.2f,%.3f,%.3f,%.3f,%.2f';

C = reshape(C,[size(X,1) nchan size(C,4)]);
C = permute(C,[1 3 2]);
C = flatten(C,2:3);

for i = 1:nchan
    chanlab1 = chanlab;
    cstr = sprintf('C%d',i);
    for j = 1:length(chanlab)
        chanlab1{j} = [chanlab1{j} cstr];
    end
    lab = [lab chanlab1];
    tplt = [tplt ',' Ctplt];
end

fid = fopen(outfile,'w');

fprintf(fid, strcat(repmat('%s,',[1 length(lab)-1]),'%s\n'),lab{:});

tplt = [tplt '\n'];
fprintf(fid, tplt, [X C]');
fclose(fid);



