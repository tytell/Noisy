function [burstnum,burstpulse] = get_burst_pulse_num(burston,burstoff, pulset, pulsedur)

burstnum = NaN(size(burston));
burstpulse = NaN(size(burston));
for i = 1:length(pulset)
    for c = 1:size(burston,2)
        indbefore = last(burstoff(:,c) < pulset(i)-pulsedur/2);
        indafter = first(burston(:,c) > pulset(i)+pulsedur/2);
        
        burstnum(indbefore+1:indafter-1,c) = 0;
        burstnum(indafter+(0:9),c) = 1:10;
        
        burstpulse(indbefore+1:indafter+9,c) = i;
    end
end
