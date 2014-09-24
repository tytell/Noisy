function data = getEntrainSpikes(data)

[spikeind,threshold] = findspikes(data.chan);

for i = 1:size(data.chan,2)
    spiket1{i} = data.t(spikeind{i});
    spikeamp1{i} = data.chan(spikeind{i},i);
end

data.spiket = catuneven(2,spiket1{:});
data.spikeamp = catuneven(2,spikeamp1{:});
data.spikethresh = threshold;

