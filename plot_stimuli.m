function plot_stimuli(data, varargin)

assert(ismember(data.StimulusType, {'Pulses','Shifts'}));

opt.phase = [];
opt.direction = [];
opt.startphase = [];
opt.phasechange = [];
opt.channel = [];
opt.raster = true;
opt.raw = true;
opt.bursts = true;
opt.meanburstphase = true;

if (isnumeric(varargin{1}))
    ind = varargin{1};
    p = 2;
else
    ind = [];
    p = 1;
end

opt = parsevarargin(opt,varargin(p:end), 2);

if (isempty(ind))
    switch data.StimulusType
        case 'Pulses'
            if isempty(opt.phase)
                opt.phase = [0 1];
            elseif numel(opt.phase) == 1
                opt.phase = [opt.phase opt.phase];
            end
            
            if isempty(opt.direction)
                opt.direction = [-1 1];
            elseif numel(opt.direction) == 1
                opt.direction = [opt.direction opt.direction];
            end
            
            ind = find((data.Phase >= opt.phase(1)) & (data.Phase <= opt.phase(2)) & ...
                (data.Direction >= opt.direction(1)) & (data.Direction <= opt.direction(2)));
            
        case 'Shifts'
            if isempty(opt.startphase) && ~isempty(opt.phase)
                opt.startphase = opt.phase;
            end
            if isempty(opt.startphase)
                opt.startphase = [0 1];
            end
            if numel(opt.startphase) == 1
                opt.startphase = opt.startphase([1 1]);
            end

            if isempty(opt.phasechange)
                opt.phasechange = [0 1];
            end
            if numel(opt.phasechange) == 1
                opt.phasechange = opt.phasechange([1 1]);
            end
            
            ind = find((data.Phase >= opt.startphase(1)) & (data.Phase <= opt.startphase(2)) & ...
                (data.PhaseChange >= opt.phasechange(1)) & (data.PhaseChange <= opt.phasechange(2)));
    end
end    

if isempty(opt.channel)
    opt.channel = 1:size(data.sig,2);
end

nchan = length(opt.channel);
nrep = length(ind);

h = -1*ones(nchan,1);

if opt.meanburstphase
    %most common maximum cycle number
    ncycle = mode(flatten(max(data.burstcyclestim)));
    
    mnphase = zeros(1,nchan);
    for c = 1:nchan
        burstcycle1 = squeeze(data.burstcyclestim(:,c,:));
        burstphase1 = squeeze(data.burstphasestim(:,c,:));
        goodcyc = (burstcycle1 < -1) | ...
            ((burstcycle1 > 4) & (burstcycle1 <= ncycle-1));
        mnphase1 = angmean(2*pi*burstphase1(goodcyc));
        mnphase(c) = mod(mnphase1/(2*pi),1);
    end
end

if size(data.tstim,3) > 1
    t1 = squeeze(data.tstim(:,1,ind));
else
    t1 = data.tstim;
end

clf;
for i = 1:nchan
    c = opt.channel(i);
    
    h(i) = subplot(nchan+1,1, i);
    sig1 = squeeze(data.sigstim(:,c,ind));
        
    if opt.raw
        good = any(isfinite(sig1),2);
        plot(t1(good,:), sig1(good,:));
    end
    
    d = max(abs(sig1(:)));
    spiket1 = squeeze(data.spiketstim(:,c,ind));
    spikey1 = d + (0:nrep-1)*d/5;

    if opt.raster
        addraster(spiket1, spikey1);
    end

    if opt.bursts
        burst1 = cat(2,data.burstonstim(:,c,ind),...
            data.burstonstim(:,c,ind)+data.burstdurstim(:,c,ind),...
            NaN(size(data.bursttstim(:,c,ind))));
        burst1 = permute(burst1,[2 1 3]);
        burst1 = flatten(burst1,1:2);
        
        bursty1 = repmat(spikey1,[size(burst1,1) 1]);
        
        addplot(burst1,bursty1,'k-', 'LineWidth',1);
        
        burstctr1 = squeeze(data.bursttstim(:,c,ind));
        bursty1 = repmat(spikey1,[size(burstctr1,1) 1]);
        addplot(burstctr1,bursty1, 'ko', 'MarkerFaceColor','w');
        if opt.meanburstphase
            axis tight;
            xl = xlim;
            stimper = 1/data.SinFreqStartHz;
            xl = round(xl/stimper);
            
            mnphase1 = (xl(1)+mnphase(i)):1:xl(2);
            if strcmp(data.StimulusType,'Shifts')
                isbefore = mnphase1 < 0;
                isduring = (mnphase1 >= 0) & (mnphase1 <= 1);
                mnphase1(isbefore) = mnphase1(isbefore) + data.PhaseChange(ind(1));
                vertplot(mnphase1(~isduring)*stimper,'k-');
                vertplot(mnphase1(isduring)*stimper,'k--');
                vertplot((mnphase1(isduring) + data.PhaseChange(ind(1)))*stimper,'k--');
            else
                mnphase1 = mnphase1 * stimper;
                vertplot(mnphase1,'k-');
            end
        end
    end
    
    axis tight;
end

h(i+1) = subplot(nchan+1,1, nchan+1);
ang1 = squeeze(data.angstim(:,1,ind));
good = any(isfinite(ang1),2);
plot(data.tstim(good), ang1(good,:));

linkaxes(h,'x');