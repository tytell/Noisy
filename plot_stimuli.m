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
opt.spacing = 0;
opt.usenewphase = true;

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
            
            if opt.usenewphase
                isphase = anginrange(2*pi*data.burstphaseatstim, 2*pi*opt.phase(1), 2*pi*opt.phase(2));
            else
                isphase = (data.Phase >= opt.phase(1)) & (data.Phase <= opt.phase(2));
            end
            ind = find(isphase & ...
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
    
    if isempty(ind)
        error('No stimuli found!');
    else
        fprintf('Found stimuli: ');
        fprintf('%d ',ind);
        fprintf('\n');
        
        [mnburstphase,~,stdburstphase] = angmean(2*pi*data.burstphaseatstim(ind));
        mnburstphase = mod(mnburstphase/(2*pi),1);
        stdburstphase = stdburstphase/(2*pi);        
        fprintf('Mean burst phase: %f +- %f\n', mnburstphase,stdburstphase);
    end
end    

if isempty(opt.channel)
    opt.channel = find(data.goodchan);
end

nchan = length(opt.channel);
nrep = length(ind);

h = -1*ones(nchan,1);

bursttype1 = data.burstprepoststim(:,opt.channel,ind);
burstper1 = data.burstperstim(:,opt.channel,ind);
burstper = nanmedian(burstper1(bursttype1 == -1));
if opt.meanburstphase
    if (data.amp > 0)
        per = 1/data.SinFreqStartHz;
    else
        per = burstper;
    end
    
    burstphase1 = data.bursttstim(:,opt.channel,ind) / per;
    burstphase1(bursttype1 ~= -1) = NaN;

    mnphase = angmean(2*pi*flatten(burstphase1,[1 3]));
    mnphase = mod(mnphase/(2*pi),1);
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
    mid1 = diff(prctile(sig1(:),[5 95]));
    if opt.raw
        good = any(isfinite(sig1),2);
        plot(repmat(t1(good,:),[1 nrep]), sig1(good,:) + ...
            repmat(0:length(ind)-1,[sum(good) 1])*mid1*opt.spacing);
    end
    
    d = max(abs(sig1(:)));
    spiket1 = squeeze(data.spiketstim(:,c,ind));
    %spiket1 = spiket1 - repmat(toff,[size(spiket1,1) 1]);
    spikey1 = mid1*opt.spacing*(length(ind)-1) + d + (0:nrep-1)*d/5;

    if opt.raster
        addraster(spiket1, spikey1);
    end
    
    if opt.bursts
        burston1 = data.burstonstim(:,c,ind);
        %burston1 = burston1 - repmat(shiftdim(toff,-1),[size(burston1,1) 1 1]);
        burstoff1 = data.burstonstim(:,c,ind)+data.burstdurstim(:,c,ind);
        %burstoff1 = burstoff1 - repmat(shiftdim(toff,-1),[size(burstoff1,1) 1 1]);
        
        burst1 = cat(2,burston1,burstoff1, NaN(size(data.bursttstim(:,c,ind))));
        burst1 = permute(burst1,[2 1 3]);
        burst1 = flatten(burst1,1:2);
        
        bursty1 = repmat(spikey1,[size(burst1,1) 1]);
        
        addplot(burst1,bursty1,'k-', 'LineWidth',1);
        
        burstctr1 = squeeze(data.bursttstim(:,c,ind));
        %burstctr1 = burstctr1 - repmat(toff,[size(burstctr1,1) 1]);
        bursty1 = repmat(spikey1,[size(burstctr1,1) 1]);
        addplot(burstctr1,bursty1, 'ko', 'MarkerFaceColor','w');
        if opt.meanburstphase
            axis tight;
            xl = xlim;
            if data.amp > 0
                stimper = 1/data.SinFreqStartHz;
            else
                stimper = burstper;
            end
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
xl = [min(data.tstim(good)) max(data.tstim(good))];

linkaxes(h,'x');
xlim(xl);