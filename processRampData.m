function processRampData(data)

nchan = size(data.spiket,2);

spikerate = zeros(length(data.Time)*3, size(data.spiket,2));
%spikerate = zeros(length(data.Time)*3, 2, size(data.spiket,2));
%spikerate = zeros(length(data.Time)*6, size(data.spiket,2));
rampvel = zeros(length(data.Time)*3);
ramppos = zeros(length(data.Time)*3);

for i = 1:length(data.Time)
    k = (i-1)*3;
    
    t1 = data.Time(i) + data.BeforeDurSec;
    
    a = first(data.t >= t1);
    b = a+1;
    while data.rampphase(b) <= 1
        b = b+1;
    end
    
    %mid = round((a+b)/2);
    
    for j = 1:nchan
        spikeind = find((data.spiket(:,j) >= data.t(a)) && ...
            (data.spiket(:,j) <= data.t(b)));
        %spikeind = find((data.spiket(:,j) >= data.t(a)) && ...
        %    (data.spiket(:,j) <= data.t(mid)));
        
        spikerate(k+1,j) = length(spikeind)/(data.t(b) - data.t(a));
    end
    rampvel(k+1) = data.RampVel(i)*data.Direction(i);
    ramppos(k+1) = data.Amplitude;
    
    a = b;
    b = a+1;
    while data.rampphase(b) <= 2
        b = b+1;
    end

    rampvel(k+2) = 0;
    
    rampvel(k+3) = -data.RampVel(i)*data.Direction(i);

end

putvar spikerate rampvel ramppos;

