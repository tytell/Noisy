function data = get_burst_params(data)

nchan = length(data.burst);
for c = 1:nchan
    on1 = data.burst(c).on;
    off1 = data.burst(c).off;
    
    int1 = zeros(size(on1));
    spikeratestd1 = zeros(size(on1));
    for i = 1:length(on1)
        isburst = (data.t >= on1(i)) & (data.t <= off1(i));
        int1(i) = trapz(data.t(isburst),abs(data.sig(isburst,c)));
        
        isburst = (data.spiket(:,c) >= on1(i)) & (data.spiket(:,c) <= off1(i));
        sr1 = 1 ./ diff(data.spiket(isburst,c));
        spikeratestd1(i) = nanstd(sr1);
    end
    
    dur1 = off1 - on1;
    spikerate1 = data.burst(c).nspike ./ dur1;
    
    data.burst(c).dur = dur1;
    data.burst(c).spikerate = spikerate1;
    data.burst(c).int = int1;
    data.burst(c).spikeratestd = spikeratestd1;
end

    
        