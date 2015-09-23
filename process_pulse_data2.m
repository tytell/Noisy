function data = process_pulse_data2(data, varargin)

opt.npre = 4;
opt.npost = 4;
opt.sampfreqlo = 50;
opt.smoothdur = 0.3;
opt.nrand = 2;
opt.showdiagnostics = true;
opt.debug = false;
opt.offsetpersegment = 0.01;
opt.savediagnostics = true;
opt.diagnosticname = '';
opt = parsevarargin(opt,varargin,2);

%get the pulse center times
pulset = data.Time + data.BeforeDurSec + data.Phase/data.stimfreq;
pulsedur2 = data.Duration(1) / 2;

nchan = size(data.sig,2);
sampfreq = 1./(data.t(2) - data.t(1));

burstper = NaN(size(data.burstt));
burstper(2:end,:) = diff(data.burstt);

%get the stimulus phase relative to the pulse
%and the burst relative phases
bt = data.burstt;
bon = data.burston;
boff = data.burstoff;

stimphase1 = NaN(length(pulset),nchan,opt.npre);
relphasemn = zeros(length(pulset),nchan);
relphasestd = zeros(length(pulset),nchan);

indclosest = zeros(length(pulset),nchan);
bursttpre = NaN(length(pulset),nchan,opt.npre);
burstonpre = NaN(length(pulset),nchan,opt.npre);
burstoffpre = NaN(length(pulset),nchan,opt.npre);
burstperpre = NaN(length(pulset),nchan,opt.npre);
bursttduring = NaN(length(pulset),nchan);
burstonduring = NaN(length(pulset),nchan);
burstoffduring = NaN(length(pulset),nchan);
burstperduring = NaN(length(pulset),nchan);
bursttpost = NaN(length(pulset),nchan,opt.npost);
burstonpost = NaN(length(pulset),nchan,opt.npost);
burstoffpost = NaN(length(pulset),nchan,opt.npost);
burstperpost = NaN(length(pulset),nchan,opt.npre);

for c = 1:nchan
    if ~data.goodchan(c)
        continue
    end

    good = ~isnan(bt(:,c));
    if all(~good)
        continue;
    end
    
    boff1 = boff(good,c);
    bon1 = bon(good,c);
    bt1 = bt(good,c);
    per1 = burstper(good,c);
        
    for i = 1:length(pulset)
        %get the burst that's closest to the pulse
        [~,indclosest(i,c)] = min(abs(bt1 - pulset(i)));
        indpre = indclosest(i,c)-1;
        
        %then get the opt.npre before
        k = indpre + (-opt.npre+1:0);

        bursttpre(i,c,:) = bt1(k);
        burstperpre(i,c,:) = per1(k);
        burstoffpre(i,c,:) = boff1(k);
        burstonpre(i,c,:) = bon1(k);

        %get the next one after the closest
        indpost = indclosest(i,c)+1;
        
        %and the opt.npost after
        k = indpost + (0:opt.npost-1);
        
        bursttpost(i,c,:) = bt1(k);
        burstonpost(i,c,:) = bon1(k);
        burstoffpost(i,c,:) = boff1(k);
        burstperpost(i,c,:) = per1(k);
        
        %and any during the burst itself
        k = indclosest(i,c);
        bursttduring(i,c) = bt1(k);
        burstonduring(i,c) = bon1(k);
        burstoffduring(i,c) = boff1(k);
        burstperduring(i,c) = per1(k);
    end
end

%most stable channel
chanstd = nanstd(flatten(burstperpre,[1 3]));
maybebad = false(size(chanstd));
for c = 1:nchan
    if ~data.goodchan(c)
        continue;
    end
    %check for any that are particularly unsteady
    if chanstd(c) > 2*nanmean(chanstd([1:c-1 c+1:end]))
        maybebad(c) = true;
    end
end

%remove any particularly unsteady channels
if any(maybebad)
    figureseries('Check freq');
    clf;
    plot(flatten(bursttpre(:,data.goodchan,:),[1 3]),...
        flatten(1./burstperpre(:,data.goodchan,:),[1 3]),'o');
    if data.amp > 0
        horizplot(data.stimfreq,'k--');
    end
    legend(data.channelnames(data.goodchan));
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    
    for c = find(maybebad)
        exc = inputyn(sprintf('%s (channel %d) seems to be highly variable. Exclude? ',data.channelnames{c},c), ...
            'default',false);
        data.goodchan(c) = ~exc;
    end
    chanstd(~data.goodchan) = NaN;
end
[~,stablechan] = min(chanstd);

%****************************************
%get the change in phase

dphasepost = NaN(length(pulset),nchan,opt.npost);
dphasepre = NaN(length(pulset),nchan,opt.npre-1);
donphasepost = NaN(length(pulset),nchan,opt.npost);
doffphasepost = NaN(length(pulset),nchan,opt.npost);
donphasepre = NaN(length(pulset),nchan,opt.npre-1);
doffphasepre = NaN(length(pulset),nchan,opt.npre-1);

for i = 1:length(pulset)
    %mean period across all channels
    perk = nanmean(flatten(burstperpre(i,:,:)));
    
    for c = 1:nchan
        %each earlier burst precedes the pulse by an integer number
        %of cycles, plus the pulse phase
        ph1 = (pulset(i) - bursttpre(i,c,:)) / perk;
        nper = round((bursttpre(i,c,end) - bursttpre(i,c,1:end-1)) / perk);
        
        %subtract off the integer number of cycles
        ph1(1:end-1) = ph1(1:end-1) - nper;
        
        %these should all be the same
        stimphase1(i,c,1:opt.npre) = ph1;

        %look for the time difference for the bursts after the stimulus
        for j = 1:opt.npost
            dt1 = bursttpost(i,c,j) - bursttpre(i,c,end) - ...
                (1+j)*perk;
            
            %average the time difference (not the phase difference), then
            %divide by the period to get the phase difference
            dphasepost(i,c,j) = nanmean(dt1)/perk;
            
            %on phase difference
            dt1 = burstonpost(i,c,j) - burstonpre(i,c,end) - ...
                (1+j)*perk;
            donphasepost(i,c,j) = nanmean(dt1)/perk;

            %off phase difference
            dt1 = burstoffpost(i,c,j) - burstoffpre(i,c,end) - ...
                (1+j)*perk;
            doffphasepost(i,c,j) = nanmean(dt1)/perk;
        end
        
        %time differences for bursts before the stimulus. Serves as a
        %control for how much the bursts naturally vary
        for j = 1:opt.npre-1
            dt1 = bursttpre(i,c,end) - bursttpre(i,c,j) - (opt.npre-j)*perk;
            dphasepre(i,c,j) = dt1 / perk;
            dt1 = burstonpre(i,c,end) - burstonpre(i,c,j) - (opt.npre-j)*perk;
            donphasepre(i,c,j) = dt1 / perk;
            dt1 = burstoffpre(i,c,end) - burstoffpre(i,c,j) - (opt.npre-j)*perk;
            doffphasepre(i,c,j) = dt1 / perk;
        end
        
        %check all the estimates
        if opt.debug && data.goodchan(c)
            clf;
            istrial = (data.t >= bursttpre(i,c,1)) & (data.t <= bursttpost(i,c,end));
            plot(data.t(istrial), data.sig(istrial,c),'k-');
            axis tight;
            overlayplot(data.t(istrial), data.ang(istrial),'r-');

            overlayplot bottom;
            yl = get(gca,'YLim');
            yy = yl(2) - 0.2*diff(yl);
            addplot(squeeze(cat(1,burstonpre(i,c,:),burstoffpre(i,c,:))),...
                yy*ones(2,opt.npre), 'b-', ...
                squeeze(cat(1,burstonduring(i,c,:),burstoffduring(i,c,:))),...
                yy*ones(2,size(burstonduring,3)), 'r-', ...
                squeeze(cat(1,burstonpost(i,c,:),burstoffpost(i,c,:))),...
                yy*ones(2,opt.npost), 'g-', 'LineWidth',2);
            
            yy = yy+0.05*diff(yl);
            text(squeeze(bursttpre(i,c,:)),yy*ones(opt.npre,1), ...
                num2str((-opt.npre:-1)'),'Color','b');
            text(squeeze(bursttpost(i,c,:)),yy*ones(opt.npost,1), ...
                num2str((1:opt.npost)'),'Color','g');
                
            vertplot(bursttpre(i,c,end)+(2:opt.npost+1)*perk,'k--');
            %text(bursttpre(i,c,end)+(nduring + (1:opt.npost))'*perk, yy*ones(opt.npost,1), ...
            %    num2str((1:opt.npost)'), 'Color','k');
            
            addplot([1;1]*(bursttpre(i,c,end)+(2:opt.npost+1)*perk) + ...
                cat(1,zeros(1,opt.npost), squeeze(dphasepost(i,c,:))'*perk), ...
                yy*ones(2,opt.npost),'k-');
            
            overlayplot top;
            pause;
        end
    end
end

%***********************************
%get change in burst duration
burstdurpre = burstoffpre - burstonpre;
burstdurduring = burstoffduring - burstonduring;
burstdurpost = burstoffpost - burstonpost;

burstdurmn = nanmean(burstdurpre,3);
burstdurfracpre = burstdurpre ./ repmat(burstdurmn,[1 1 opt.npre]);
burstdurfracduring = burstdurduring ./ burstdurmn;
burstdurfracpost = burstdurpost ./ repmat(burstdurmn,[1 1 opt.npost]);

%***********************************
%get the relative phases of the burst channels
%relative to the most stable channel
for i = 1:length(pulset)
    c = stablechan;
    bt1 = bursttpre(i,c,:);
    for c2 = 1:nchan
        if (c2 == c) || ~data.goodchan(c)
            continue;
        end
        bt2 = bursttpre(i,c2,:);
        
        %matrix of all the differences
        d = repmat(bt2(:)',[length(bt1) 1]) - repmat(bt1(:),[1 length(bt2)]);
        phdiff1 = d/perk;
        
        [phdiffmn1,~,phdiffstd1] = angmean(2*pi*phdiff1(:));
        relphasemn(i,c2) = mod(phdiffmn1/(2*pi),1);
        relphasestd(i,c2) = phdiffstd1/(2*pi);
    end
end

%stimphase1 is the stim phase relative to bursts for each channel
%add on the relative phase between the channels to get a consistent phase
%relative to the most stable channel
stimphase = stimphase1 + repmat(relphasemn,[1 1 opt.npre]);
[stimphase,~,stimphasestd] = angmean(2*pi*flatten(stimphase,2:3), 2);
stimphase = stimphase/(2*pi);
stimphasestd = stimphasestd/(2*pi);

%now adjust the stimphase so that it's relative to the left side burst at
%the stimulus position
tok = regexp(data.channelnames{stablechan},'(L|R)(\d+)','once','tokens');
if isempty(tok)
    error('Could not parse channel name %s',data.channelnames{stablechan});
end
if tok{1} == 'R'
    stimphase = stimphase + 0.5;
end
pos = str2double(tok{2});
dist = data.stimuluspos - pos;
stimphase = stimphase - dist * opt.offsetpersegment;
stimphase = mod(stimphase,1);

if isfield(data,'stimcyclet');
    data = rmfield(data,{'stimcyclet','burstspercycle','burstcycle','burstcyclet'});
end

%save the values in the data structure
data.stimphase = stimphase;
data.stimphasestd = stimphasestd;
data.stablechan = stablechan;
data.pulset = pulset;
data.burstfreq = 1./burstper;
data.bursttpre = bursttpre;
data.burstonpre = burstonpre;
data.burstoffpre = burstoffpre;
data.burstperpre = burstperpre;
data.bursttduring = bursttduring;
data.burstonduring = burstonduring;
data.burstoffduring = burstoffduring;
data.burstperduring = burstperduring;
data.bursttpost = bursttpost;
data.burstonpost = burstonpost;
data.burstoffpost = burstoffpost;
data.burstperpost = burstperpost;
data.dphasepost = dphasepost;
data.dphasepre = dphasepre;
data.donphasepost = donphasepost;
data.doffphasepost = doffphasepost;
data.donphasepre = donphasepre;
data.doffphasepre = doffphasepre;
data.burstdurfracpre = burstdurfracpre;
data.burstdurfracduring = burstdurfracduring;
data.burstdurfracpost = burstdurfracpost;

if opt.showdiagnostics 
    %histogram of stimulus phases
    fig(1) = figureseries('Stimulus phase');
    clf;
    edges = 0:1/6:1;
    isleft = data.Direction == -1;
    nleft = histc(stimphase(isleft),edges);
    nright = histc(stimphase(~isleft),edges);
    ctrs = (edges(1:end-1) + edges(2:end))/2;
    bar(ctrs,[nleft(1:end-1) nright(1:end-1)],1,'stacked');
    xlabel('Phase');
    ylabel('Number of pulses');
    legend('left','right');
    input('Hit return to continue.');
    
    [bursttnear,burstfreqnear] = eventtrigger(data.pulset, data.burstt', [-2 8], ...
        data.burstfreq', 'eventrange');
    bursttnear = bursttnear - repmat(data.pulset',[size(bursttnear,1) 1 nchan]);
    
    %burst frequency before and after the stimulus
    fig(2) = figureseries('Frequency');
    clf;
    channear = repmat(shiftdim(1:nchan,-1),[size(bursttnear,1) size(bursttnear,2) 1]);
    burstnumnear = repmat((-2:8)',[1 size(bursttnear,2) nchan]);
    
    h = plotgroups(burstnumnear(:),burstfreqnear(:),{channear(:)}, {'cmf'},'means','xoff',0);
    if data.amp > 0
        horizplot(data.stimfreq,'k--');
    end
    xlabel('Burst number');
    ylabel('Frequency (Hz)');
    legend(h(data.goodchan),data.channelnames{data.goodchan});
    input('Hit return to continue.');
    
    %raster plot of all the bursts and stimuli
    spiket1 = eventtrigger(data.pulset, data.spiket', [-4 5],'timerange');
    
    [~,ord] = sortrows([data.Direction data.stimphase]);

    % t0 = burstoff1(:,:,1);
    % t0 = t0(burstnum1(:,:,1) == -1);
    if data.amp > 0
        basefreq = repmat(data.stimfreq, size(spiket1));
    else
        basefreq = nanmean(flatten(1./burstperpre,2:3),2);
        basefreq = repmat(shiftdim(basefreq,-1), [size(spiket1,1) 1 nchan]);
    end
    t0 = bursttpre(:,stablechan,end); % + relphasemn./squeeze(basefreq(1,:,:));
    %t0 = nanmean(t0,2);
    %t0 = data.pulset - stimphase./basefreq(1,:,1)';
    
    spikeph1 = (spiket1 - repmat(t0', [size(spiket1,1) 1 nchan])) .* basefreq;
    spikeph1 = spikeph1(:,ord,:);
    
    stimphctr1 = (data.pulset - t0).*basefreq(1,:,1)';
    stimphend1 = (data.pulset + pulsedur2 - t0).*basefreq(1,:,1)';
    
    burstphpre1 = (bursttpre - repmat(t0, [1 nchan opt.npre])) .* ...
        repmat(squeeze(basefreq(1,:,:)), [1 1 opt.npre]);
    
    fig(3) = figureseries('Raster');
    clf;
    
    gc = find(data.goodchan);
    clear h;
    for i = 1:length(gc)
        c = gc(i);
        subplot(length(gc),1,i);
        
        raster(spikeph1(:,:,c),'k');
        
        h(1) = addplot(stimphctr1(ord),1:length(ord),'r*-');
        addplot(stimphend1(ord),1:length(ord),'r-');
        
        h(2) = addplot(burstphpre1(ord,c,end), 1:length(ord), 'b-');
        
        h(3) = addplot(burstphpre1(ord,c,end) + 2, 1:length(ord),'b--');
        h(4) = addplot(burstphpre1(ord,c,end) + 2 + dphasepost(ord,c,1), ...
            1:length(ord),'g-');

        xlabel('Phase');
        ylabel('Trial');
        if i == 1
            legend(h,'stim','burst -1','predicted burst 1','actual burst 1');
        end
    end
    input('Hit return to continue.');
    
    %phase response curves for the first two bursts after the stimulus
    fig(4) = figureseries('PRC');
    hax(1) = subplot(2,2,1);
    isleft = data.Direction == -1;
    h1 = plot(stimphase(isleft),dphasepre(isleft,:,end),'b.');
    h2 = addplot(stimphase(isleft),dphasepost(isleft,:,1),'ko');
    legend([h1(1) h2(1)],'No pulse','With pulse');
    xlabel('Phase');
    ylabel({'Burst 1','Change in phase'});
    title('Left pulses');
    yl(1,:) = get(gca,'YLim');
    
    hax(2) = subplot(2,2,2);
    plot(stimphase(~isleft),dphasepre(~isleft,:,end),'b.');
    addplot(stimphase(~isleft),dphasepost(~isleft,:,1),'ko');
    xlabel('Phase');
    ylabel({'Burst 1','Change in phase'});
    title('Right pulses');
    yl(2,:) = get(gca,'YLim');

    hax(3) = subplot(2,2,3);
    plot(stimphase(isleft),dphasepre(isleft,:,end-1),'b.');
    addplot(stimphase(isleft),dphasepost(isleft,:,2),'ko');
    xlabel('Phase');
    ylabel({'Burst 2','Change in phase'});
    yl(3,:) = get(gca,'YLim');
    
    hax(4) = subplot(2,2,4);
    plot(stimphase(~isleft),dphasepre(~isleft,:,end-1),'b.');
    addplot(stimphase(~isleft),dphasepost(~isleft,:,2),'ko');
    xlabel('Phase');
    ylabel({'Burst 2','Change in phase'});
    yl(4,:) = get(gca,'YLim');
    
    set(hax,'YLim',max(yl));
    input('Hit return to continue.');
    
    %burst durations
    fig(5) = figureseries('Duration');
    hax(1) = subplot(2,2,1);
    isleft = data.Direction == -1;
    h1 = plot(stimphase(isleft),burstdurfracpre(isleft,:,end),'b.');
    h2 = addplot(stimphase(isleft),burstdurfracduring(isleft,:),'ko');
    legend([h1(1) h2(1)],'No pulse','With pulse');
    xlabel('Phase');
    ylabel({'During','Relative burst duration'});
    title('Left pulses');
    yl(1,:) = get(gca,'YLim');
    
    hax(2) = subplot(2,2,2);
    plot(stimphase(~isleft),burstdurfracpre(~isleft,:,end),'b.');
    addplot(stimphase(~isleft),burstdurfracduring(~isleft,:),'ko');
    xlabel('Phase');
    ylabel({'During','Relative burst duration'});
    title('Right pulses');
    yl(2,:) = get(gca,'YLim');

    hax(3) = subplot(2,2,3);
    plot(stimphase(isleft),burstdurfracpre(isleft,:,end-1),'b.');
    addplot(stimphase(isleft),burstdurfracpost(isleft,:,1),'ko');
    xlabel('Phase');
    ylabel({'Burst 1','Relative burst duration'});
    yl(3,:) = get(gca,'YLim');

    hax(4) = subplot(2,2,4);
    plot(stimphase(~isleft),burstdurfracpre(~isleft,:,end-1),'b.');
    addplot(stimphase(~isleft),burstdurfracpost(~isleft,:,1),'ko');
    xlabel('Phase');
    ylabel({'Burst 1','Relative burst duration'});
    yl(4,:) = get(gca,'YLim');
    
    set(hax,'YLim',max(yl));
    input('Hit return to continue.');
    
    if opt.savediagnostics
        fprintf('Saving figures...\n');
        for i = 1:5
            nm = sprintf('%s-%d%s',opt.diagnosticname(1:end-4),i,opt.diagnosticname(end-3:end));
            
            fprintf('%d...\n',i);
            print(fig(i),'-dpdf',nm);
        end
        fprintf('Done\n');
    end
end


