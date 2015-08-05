function check_data(data)

cycrange = [-4 6];
trange = cycrange / data.stimfreq;
dt = data.t(2) - data.t(1);

col = 'krgbm';

nchan = size(data.burstt,2);
nstim = length(data.Time);

bursttype = NaN(size(data.burstt));
burstindpre = zeros(4,nstim,nchan);
burstindpost = zeros(6,nstim,nchan);

stimt = (trange(1):dt:trange(2))';
stimt = repmat(stimt,[1 nstim]);
stim = NaN(size(stimt,1),nstim);

indgoodchan = find(data.goodchan);

stimctr = data.Time + data.Phase/data.stimfreq + data.BeforeDurSec;

spiket = cell(nstim,nchan);
for i = 1:length(data.Time)
    for c = 1:nchan
        if data.goodchan(c)
            indpre = find(data.burstoff(:,c) < stimctr(i) - data.Duration(i)/2, ...
                4,'last');
            bursttype(indpre) = -4:-1;
            burstindpre(:,i,c) = indpre;
            
            inddur = find((data.burstoff(:,c) >= stimctr(i) - data.Duration(i)/2) & ...
                (data.burston(:,c) <= stimctr(i) + data.Duration(i)/2));
            bursttype(inddur) = 0;
            
            indpost = find((data.burston(:,c) > stimctr(i) + data.Duration(i)/2), ...
                6, 'first');
            bursttype(indpost) = 1:6;
            burstindpost(:,i,c) = indpost;
            
            isspike = (data.spiket(:,c) >= data.burston(indpre(1),c)) & ...
                (data.spiket(:,c) <= data.burstoff(indpost(end),c));
            spiket{i,c} = data.spiket(isspike,c);
        else
            spiket{i,c} = NaN;
        end
    end
    stimt(:,i) = stimt(:,i) + stimctr(i);
    stim(:,i) = interp1(data.t,data.ang, stimt(:,i));
end

spiket = catuneven(3,spiket{:});
spiket = reshape(spiket,[],nstim,nchan);
nspike = size(spiket,1);

burstt0 = data.burstt(burstindpre(end,:,2),2)';
t0 = data.Time' + data.BeforeDurSec';
%t0 = burstt0;

stimtoff = data.Time' + data.BeforeDurSec + data.Phase'/data.stimfreq - burstt0;
stimphase = stimtoff*data.stimfreq;

spikeph = (spiket - repmat(t0,[nspike 1 nchan])) * data.stimfreq;

burstphmn = NaN(1,nchan);
for c = 1:nchan
    if data.goodchan(c)
        burstphmn(c) = angmean(2*pi*data.burstphase(isnan(bursttype(:,c)),c)) / (2*pi);
    end
end

[~,ord] = sortrows([data.Direction,stimphase']);

figure(1);
clf;

h(1) = subplot(4,1,1:3);
hold on;
for c = 1:nchan
    if data.goodchan(c)
        raster(spikeph(:,ord,c)-burstphmn(c), repmat(1:length(ord), [nspike 1]), col(c));
    end
end

c = first(data.goodchan == 1);
addplot((data.Time(ord)' + data.BeforeDurSec - t0(ord)) * data.stimfreq + data.Phase(ord)' - burstphmn(c), 1:length(ord), 'k*');

h(2) = subplot(4,1,4);
plot(stimt - repmat(t0,[size(stimt,1) 1]), stim);

linkaxes(h,'x');

jitter = linspace(-0.03,0.03,4)';
figure(3);
clf;

h(1) = subplot(2,1,1);
hold(h(1),'on');
h(2) = subplot(2,1,2);
hold(h(2),'on');

burstonphase = NaN(size(data.burston));
good = isfinite(data.burston);
burstonphase(good) = interp1(data.t,data.phase, data.burston(good));

burstphpre = NaN(size(burstindpre));
burstphpost = NaN(size(burstindpost));
for c = 1:nchan
    if data.goodchan(c)
        burstphpre1 = burstonphase(burstindpre(:,:,c),c) - burstphmn(c);
        burstphpre1 = mod(burstphpre1 + 0.5, 1) - 0.5;
        burstphpre(:,:,c) = reshape(burstphpre1,[size(burstphpre,1) nstim]);
        burstphpost1 = burstonphase(burstindpost(:,:,c),c) - burstphmn(c);
        burstphpost1 = mod(burstphpost1 + 0.5, 1) - 0.5;
        burstphpost(:,:,c) = reshape(burstphpost1,[size(burstphpost,1) nstim]);

        dirs = [-1 1];
        for d = 1:2
            isdir = data.Direction == dirs(d);
            plot(h(d), repmat(data.Phase(isdir)',[4 1]) + repmat(jitter(1:4), [1 sum(isdir)]), ...
                [burstphpre(end,isdir,c); burstphpost(1:3,isdir,c)], '-', 'Color',[0.7 0.7 0.7]);
            plot(h(d), data.Phase(isdir)+jitter(1), burstphpre(:,isdir,c), 'o', 'Color',[0.7 0.7 0.7]);
            plot(h(d), data.Phase(isdir)+jitter(2), burstphpost(1,isdir,c), 'r*');
            plot(h(d), data.Phase(isdir)+jitter(3), burstphpost(2,isdir,c), 'gd');
            plot(h(d), data.Phase(isdir)+jitter(4), burstphpost(3,isdir,c), 'bs');
        end
    end
end
hold off;

%addplot(data.Phase+jitter(3), burstphpost(2,:), 'g+');
%addplot(data.Phase+jitter(4), burstphpost(3:end,:), 'bd');

% figure(4);
% plot(stimphase, burstphpre, 'o', 'Color',[0.7 0.7 0.7]);
% addplot(stimphase, burstphpost(1,:), 'r*');
% addplot(stimphase, burstphpost(2,:), 'g+');
% addplot(stimphase, burstphpost(3:end,:), 'bd');
% 
% 




