function plot_stimuli(data, varargin)

assert(ismember(data.StimulusType, {'Pulses'}));

opt.phase = [];
opt.direction = [];
opt.channel = [];
opt.raster = true;
opt.raw = true;
opt.bursts = true;
opt.meanburstphase = true;
opt.spacing = 0;

if (isnumeric(varargin{1}))
    ind = varargin{1};
    p = 2;
else
    ind = [];
    p = 1;
end

opt = parsevarargin(opt,varargin(p:end), 2);

if (isempty(ind))
    if isempty(opt.phase)
        opt.phase = [0 1];
    elseif numel(opt.phase == 1)
        opt.phase = [opt.phase opt.phase];
    end

    if isempty(opt.direction)
        opt.direction = [-1 1];
    elseif numel(opt.direction == 1)
        opt.direction = [opt.direction opt.direction];
    end
    
    ind = find((data.Phase >= opt.phase(1)) & (data.Phase <= opt.phase(2)) & ...
        (data.Direction >= opt.direction(1)) & (data.Direction <= opt.direction(2)));
    if isempty(ind)
        error('No stimuli found!');
    else
        fprintf('Found stimuli: ');
        fprintf('%d ',ind);
        fprintf('\n');
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

clf;
for i = 1:nchan
    c = opt.channel(i);

    
    h(i) = subplot(nchan+1,1, i);
    sig1 = squeeze(data.sigstim(:,c,ind));
    mid1 = diff(prctile(sig1(:),[5 95]));
    if opt.raw
        good = any(isfinite(sig1),2);
        plot(data.tstim(good), sig1(good,:) + ...
            repmat(0:length(ind)-1,[sum(good) 1])*mid1*opt.spacing);
    end
    
    d = max(abs(sig1(:)));
    spiket1 = squeeze(data.spiketstim(:,c,ind));
    spikey1 = mid1*opt.spacing*(length(ind)-1) + d + (0:nrep-1)*d/5;

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
            xl = round(xl/stimper)*stimper;
            
            mnphase1 = (xl(1)+mnphase(i)*stimper):stimper:xl(2);
            vertplot(mnphase1,'k-');
        end
    end
    
    axis tight;
end

h(i+1) = subplot(nchan+1,1, nchan+1);
ang1 = squeeze(data.angstim(:,1,ind));
good = any(isfinite(ang1),2);
plot(data.tstim(good), ang1(good,:));

linkaxes(h,'x');