function processRampData(data)

dt = (data.t(2) - data.t(1));
looprate = 20;
offset = (1/looprate)/2;

nchan = size(data.spiket,2);

spikerate = zeros(length(data.Time)*3, size(data.spiket,2));
%spikerate = zeros(length(data.Time)*3, 2, size(data.spiket,2));
%spikerate = zeros(length(data.Time)*6, size(data.spiket,2));
spikenum = zeros(length(data.Time)*3, size(data.spiket, 2));
rampvel = zeros(length(data.Time)*3, 1);
ramppos = zeros(length(data.Time)*3, 1);
holdtime = zeros(length(data.Time)*3, 1);

zerospikerate = zeros(length(data.Time), size(data.spiket, 2));
zerospikenum = zeros(length(data.Time), size(data.spiket, 2));
zerorampvel = zeros(length(data.Time), 1);
zeroramppos = zeros(length(data.Time), 1);
zeroholdtime = zeros(length(data.Time), 1);


for i = 1:length(data.Time) %Steps through each start stimulus start t
    k = (i-1)*3; % Divides into the 3 phases
    
    %Phase 1 of ramp
    
    t1 = data.Time(i) + data.BeforeDurSec + offset; %data.Time doesn't include time before 1st stimulus
    t2 = data.Amplitude(i) / data.RampVel(i) + t1;
    %t2 = data.t(first(data.ang >= data.Amplitude(i)));
    
    rampdur = t2-t1;
    %mid = round((t1+t2)/2);
    
    for j = 1:nchan
        spikeind = find((data.spiket(:,j) >= t1) & (data.spiket(:,j) <= t2));
        %spikeind = find((data.spiket(:,j) >= data.t(a)) & ...
        %    (data.spiket(:,j) <= data.t(mid)));
        spikenum(k+1,j) = length(spikeind);
        spikerate(k+1,j) = length(spikeind)/(t2 - t1);
    end
    
    rampvel(k+1) = data.RampVel(i)*data.Direction(i);
    ramppos(k+1) = data.Amplitude(i)*data.Direction(i);
    holdtime(k+1) = data.Duration(i);
    
    %Phase 2 of ramp
    t1 = t2;
    t2 = t2 + data.Duration(i);
    
    for j = 1:nchan
        spikeind = find((data.spiket(:,j) >= t1) & ...
            (data.spiket(:,j) <= t2));
        spikenum(k+2,j) = length(spikeind);
        spikerate(k+2,j) = length(spikeind)/(t2 - t1);
    end
    
    rampvel(k+2) = 0;
    ramppos(k+2) = data.Amplitude(i)*data.Direction(i);
    holdtime(k+2) = data.Duration(i);
    
    %Phase 3 of ramp
    t1 = t2;
    t2 = t2 + rampdur;
    
    for j = 1:nchan
        spikeind = find((data.spiket(:,j) >= t1) & ...
            (data.spiket(:,j) <= t2));
        spikenum(k+3,j) = length(spikeind);
        spikerate(k+3,j) = length(spikeind)/(t2 - t1);
    end
    
    rampvel(k+3) = -data.RampVel(i)*data.Direction(i);
    ramppos(k+3) = data.Amplitude(i)*data.Direction(i);
    holdtime(k+3) = data.Duration(i);

    %Zero spike rate (spike rate at neutral position)
    if i == length(data.Time)
        t1 = t2;
        t2 = data.t(length(data.t));
        
         for j = 1:nchan
             spikeind = find((data.spiket(:,j) >= t1) & ...
                (data.spiket(:,j) <= t2));
             zerospikenum(i,j) = length(spikeind);
             zerospikerate(i,j) = length(spikeind)/(t2 - t1);
         end

         zerorampvel(i) = data.RampVel(i)*data.Direction(i);
         zeroramppos (i) = data.Amplitude(i)*data.Direction(i);
         zeroholdtime(i) = t2 - t1;
         
        break;
    end
    
    t1 = t2;
    t2 = data.Time(i+1) + data.BeforeDurSec + offset;

    for j = 1:nchan
        spikeind = find((data.spiket(:,j) >= t1) & ...
            (data.spiket(:,j) <= t2));
        zerospikenum(i,j) = length(spikeind);
        zerospikerate(i,j) = length(spikeind)/(t2 - t1);
    end

    zerorampvel(i) = data.RampVel(i)*data.Direction(i);
    zeroramppos (i) = data.Amplitude(i)*data.Direction(i);
    zeroholdtime(i) = t2 - t1;
    
end

putvar spikerate rampvel ramppos spikenum holdtime zerospikenum ...
    zerospikerate zerorampvel zeroramppos zeroholdtime;

