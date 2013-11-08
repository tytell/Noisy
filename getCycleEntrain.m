function data = getCycleEntrain(data, varargin)

opt.navg = 11;
opt = parsevarargin(opt,varargin, 2);

ncycle = size(data.stimcyclet,1);
spikecycle = NaN(size(data.spiket));
Rcycle = NaN(ncycle,size(data.spiket,2));
phasecycle = NaN(ncycle,size(data.spiket,2));
npercycle = NaN(ncycle,size(data.spiket,2));

navg2 = floor(opt.navg/2);
navg = opt.navg;

goodcycle = isfinite(data.cycle);
for i = 1:size(data.spiket,2)
    good = isfinite(data.spiket(:,i));
    spikecycle(good,i) = interp1(data.t(goodcycle), data.cycle(goodcycle), data.spiket(good,i));
    
    good = isfinite(spikecycle(:,i));
    c1 = accumarray(spikecycle(good,i), data.spikephase(good,i), [], ...
        @(x) sum(cos(2*pi*x)));
    s1 = accumarray(spikecycle(good,i), data.spikephase(good,i), [], ...
        @(x) sum(sin(2*pi*x)));
    n1 = accumarray(spikecycle(good,i), data.spikephase(good,i), [], ...
        @(x) sum(isfinite(x)));
    
    cblock1 = cumsum(c1);
    sblock1 = cumsum(s1);
    nblock1 = cumsum(n1);
    
    %get means in blocks of opt.navg cycles
    c2 = NaN(size(c1));
    c2(navg2+1:end-navg2-1) = cblock1(navg+1:end) - cblock1(1:end-navg);
    s2 = NaN(size(s1));
    s2(navg2+1:end-navg2-1) = sblock1(navg+1:end) - sblock1(1:end-navg);
    n2 = zeros(size(n1));
    n2(navg2+1:end-navg2-1) = nblock1(navg+1:end) - nblock1(1:end-navg);
    
    Rcycle(1:length(c1),i) = sqrt(c2.^2 + s2.^2)./n2;
    phasecycle(1:length(c1),i) = mod(atan2(s2./n2,c2./n2)/(2*pi),1);
    npercycle(1:length(c1),i) = n2 / opt.navg;
end

data.spikeRcycle = Rcycle;
data.spikephasecycle = phasecycle;
data.nspikespercycle = npercycle;

if (length(data.amp) == length(data.t))
    data.ampcycle = interp1(data.t,data.amp, data.stimcyclet);
end
if (length(data.stimfreq) == length(data.t))
    data.stimfreqcycle = interp1(data.t,data.stimfreq, data.stimcyclet);
end
if (length(data.noise) == length(data.t))
    data.noisecycle = interp1(data.t,data.noise, data.stimcyclet);
end
