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

burstfreq = NaN(size(data.bursttstim));
burstfreq(1:end-1,:,:) = 1 ./ diff(data.bursttstim,[],1);

burstint0 = catuneven(1,data.burst.int)';
burstspikerate0 = catuneven(1,data.burst.spikerate)';

burstspikerate = NaN(size(burstfreq));
good = isfinite(data.burstindstim);
burstspikerate(good) = burstspikerate0(data.burstindstim(good));

nchan = size(data.sig,2);
mnphase = zeros(1,nchan);
stdphase = zeros(1,nchan);
mnonphase = zeros(1,nchan);
stdonphase = zeros(1,nchan);
mndur = zeros(1,nchan);
stddur = zeros(1,nchan);
mnspikerate = zeros(1,nchan);
stdspikerate = zeros(1,nchan);
for c = 1:nchan
    prepost1 = squeeze(data.burstprepoststim(:,c,:));
    burstphase1 = squeeze(data.burstphasestim(:,c,:));
    goodcyc = (prepost1 == -1) | (prepost1 == 3);
    [mnphase1,~,stdphase1] = angmean(2*pi*burstphase1(goodcyc));
    mnphase(c) = mod(mnphase1/(2*pi),1);
    stdphase(c) = stdphase1/(2*pi);
    
    onphase1 = onphase(:,c,:);
    [mnonphase1,~,stdonphase1] = angmean(2*pi*onphase1(goodcyc));
    mnonphase(c) = mod(mnonphase1/(2*pi),1);
    stdonphase(c) = stdonphase1/(2*pi);
    
    dur1 = data.burstdurstim(:,c,:);
    mndur(c) = nanmean(dur1(goodcyc));
    stddur(c) = nanstd(dur1(goodcyc));
    
    spikerate1 = burstspikerate(:,c,:);
    mnspikerate(c) = nanmean(spikerate1(goodcyc));
    stdspikerate(c) = nanstd(spikerate1(goodcyc));
end

cycle = data.burstcyclestim;
phase = data.burstphasestim;
prepost = data.burstprepoststim;

expectedt = (cycle + repmat(mnphase,[size(cycle,1) 1 size(cycle,3)]))*stimper;
dphase = (data.bursttstim - expectedt)/stimper;
dphasez = dphase ./ repmat(stdphase,[size(cycle,1) 1 size(cycle,3)]);

expectedon = (cycle + repmat(mnonphase,[size(cycle,1) 1 size(cycle,3)]))*stimper;
donphase = (data.burstonstim - expectedon)/stimper;
donphasez = donphase ./ repmat(stdonphase,[size(cycle,1) 1 size(cycle,3)]);

dur = data.burstdurstim / stimper;
ddur = dur - repmat(mndur/stimper,[size(cycle,1) 1 size(cycle,3)]);
ddurz = ddur ./ repmat(stddur/stimper,[size(cycle,1) 1 size(cycle,3)]);

dburstspikerate = burstspikerate - repmat(mnspikerate,[size(cycle,1) 1 size(cycle,3)]);
dburstspikeratez = dburstspikerate ./ repmat(stdspikerate,[size(cycle,1) 1 size(cycle,3)]);

prepost = permute(prepost,[1 3 2]);
cycle = permute(cycle,[1 3 2]);
phase = permute(phase,[1 3 2]);
burstfreq = permute(burstfreq,[1 3 2]);
dphase = permute(dphase,[1 3 2]);
dphasez = permute(dphasez,[1 3 2]);
donphase = permute(donphase,[1 3 2]);
donphasez = permute(donphasez,[1 3 2]);
dur = permute(dur,[1 3 2]);
ddur = permute(ddur,[1 3 2]);
ddurz = permute(ddurz,[1 3 2]);
burstspikerate = permute(burstspikerate, [1 3 2]);
dburstspikerate = permute(dburstspikerate, [1 3 2]);
dburstspikeratez = permute(dburstspikeratez, [1 3 2]);

good = ~isnan(cycle) & ~isnan(prepost);
cpind = NaN(size(cycle));
[cp,~,cpind(good)] = unique([cycle(good) prepost(good)],'rows');
phase2 = NaN(length(cp),size(phase,2),nchan);
burstfreq2 = NaN(size(phase2));
dphase2 = NaN(size(phase2));
dphasez2 = NaN(size(phase2));
donphase2 = NaN(size(phase2));
donphasez2 = NaN(size(phase2));
dur2 = NaN(size(phase2));
ddur2 = NaN(size(phase2));
ddurz2 = NaN(size(phase2));
burstspikerate2 = NaN(size(phase2));
dburstspikerate2 = NaN(size(phase2));
dburstspikeratez2 = NaN(size(phase2));
for c = 1:nchan
    for i = 1:size(phase,2)
        for j = 1:size(cp,1)
            iscycle = cpind(:,i,c) == j;
            if sum(iscycle) == 1
                phase2(j,i,c) = phase(iscycle,i,c);
                burstfreq2(j,i,c) = burstfreq(iscycle,i,c);
                dphase2(j,i,c) = dphase(iscycle,i,c);
                dphasez2(j,i,c) = dphasez(iscycle,i,c);
                donphase2(j,i,c) = donphase(iscycle,i,c);
                donphasez2(j,i,c) = donphasez(iscycle,i,c);
                dur2(j,i,c) = dur(iscycle,i,c);
                ddur2(j,i,c) = ddur(iscycle,i,c);
                ddurz2(j,i,c) = ddurz(iscycle,i,c);
                burstspikerate2(j,i,c) = burstspikerate(iscycle,i,c);
                dburstspikerate2(j,i,c) = dburstspikerate(iscycle,i,c);
                dburstspikeratez2(j,i,c) = dburstspikeratez(iscycle,i,c);
            end
        end
    end
end

cycle2 = repmat(cp(:,1),[1 size(phase,2)]);
prepost2 = repmat(cp(:,2),[1 size(phase,2)]);

C = cat(4,phase2,burstfreq2,dphase2,dphasez2, donphase2,donphasez2, ...
    dur2,ddur2,ddurz2, burstspikerate2,dburstspikerate2,dburstspikeratez2);

stimphase = repmat(data.Phase',[size(cycle2,1) 1]);
stimdirec = repmat(data.Direction',[size(cycle2,1) 1]);
stimamp = repmat(data.Amplitude',[size(cycle2,1) 1]);
stimdur = repmat(data.Duration',[size(cycle2,1) 1]);
stimnum = repmat(1:length(data.Phase),[size(cycle2,1) 1]);

X = [stimnum(:) stimphase(:) stimdirec(:) stimamp(:) stimdur(:) cycle2(:) prepost2(:)];
lab = {'StimNum','StimPhase','StimDirec','StimAmp','StimDur','CycleNum','CycleType'};
chanlab = {'Phase','Freq','DPhase','DPhaseZ','DOnPhase','DOnPhaseZ',...
    'Dur','DDur','DDurZ','SpikeRate','DSpikeRate','DSpikeRateZ'};


tplt = '%d,%.3f,%d,%.1f,%.2f,%.1f,%d';
Ctplt = '%.3f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.3f,%.3f';

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



