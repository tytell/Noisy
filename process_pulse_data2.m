function data = process_pulse_data2(data, varargin)

opt.npre = 4;
opt.npost = 4;
opt.sampfreqlo = 50;
opt.smoothdur = 0.3;
opt.nrand = 2;
opt.showdiagnostics = true;
opt.debug = false;
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

bursttpre = NaN(length(pulset),nchan,opt.npre);
burstonpre = NaN(length(pulset),nchan,opt.npre);
burstoffpre = NaN(length(pulset),nchan,opt.npre);
burstperpre = NaN(length(pulset),nchan,opt.npre);
bursttduring = NaN(length(pulset),nchan,2);
burstonduring = NaN(length(pulset),nchan,2);
burstoffduring = NaN(length(pulset),nchan,2);
bursttpost = NaN(length(pulset),nchan,opt.npost);
burstonpost = NaN(length(pulset),nchan,opt.npost);
burstoffpost = NaN(length(pulset),nchan,opt.npost);
dphasepost = NaN(length(pulset),nchan,opt.npost);

for c = 1:nchan
    if ~data.goodchan(c)
        continue
    end

    good = ~isnan(bt(:,c));
    boff1 = boff(good,c);
    bon1 = bon(good,c);
    bt1 = bt(good,c);
    per1 = burstper(good,c);
        
    for i = 1:length(pulset)
        %get the last burst that ends before the pulse begins
        indpre = last(boff1 < pulset(i)-pulsedur2);
        
        %then get the opt.npre before
        k = indpre + (-opt.npre+1:0);

        bursttpre(i,c,:) = bt1(k);
        burstperpre(i,c,:) = per1(k);
        burstoffpre(i,c,:) = boff1(k);
        burstonpre(i,c,:) = bon1(k);

        %get the last burst that ends before the pulse begins
        indpost = first(bon1 > pulset(i)+pulsedur2);
        
        %and the opt.npost after
        k = indpost + (0:opt.npost-1);
        
        bursttpost(i,c,:) = bt1(k);
        burstonpost(i,c,:) = bon1(k);
        burstoffpost(i,c,:) = boff1(k);
        
        %and any during the burst itself
        k = indpre+1:indpost-1;
        bursttduring(i,c,1:length(k)) = bt1(k);
        burstonduring(i,c,1:length(k)) = bon1(k);
        burstoffduring(i,c,1:length(k)) = boff1(k);
    end
end
    
nduring = zeros(length(pulset),nchan);
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
        
        stimphase1(i,c,1:opt.npre) = ph1;

        nduring(i,c) = sum(isfinite(bursttduring(i,c,:)));
        nduring1 = nduring(i,c);
        for j = 1:opt.npost
            dt1 = bursttpost(i,c,j) - bursttpre(i,c,:);
            nper = nduring1 + j + (opt.npre-1:-1:0);
            
            %average the time difference (not the phase difference)
            dt1 = dt1 - shiftdim(nper,-1)*perk;
            dphasepost(i,c,j) = nanmean(dt1)/perk;
        end
        
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
                
            vertplot(bursttpre(i,c,end)+(1:nduring1+opt.npost)*perk,'k--');
            %text(bursttpre(i,c,end)+(nduring + (1:opt.npost))'*perk, yy*ones(opt.npost,1), ...
            %    num2str((1:opt.npost)'), 'Color','k');
            
            addplot([1;1]*(bursttpre(i,c,end)+(nduring1 + (1:opt.npost))*perk) + ...
                cat(1,zeros(1,opt.npost), squeeze(dphasepost(i,c,:))'*perk), ...
                yy*ones(2,opt.npost),'k-');
            
            overlayplot top;
            pause;
        end
    end
    
    c = first(data.goodchan);
    bt1 = bursttpre(i,c,:);
    for c2 = c+1:nchan
        bt2 = bursttpre(i,c2,:);
        
        %matrix of all the differences
        d = repmat(bt2(:)',[length(bt1) 1]) - repmat(bt1(:),[1 length(bt2)]);
        phdiff1 = d/perk;
        
        [phdiffmn1,~,phdiffstd1] = angmean(2*pi*phdiff1(:));
        relphasemn(i,c2) = mod(phdiffmn1/(2*pi),1);
        relphasestd(i,c2) = phdiffstd1/(2*pi);
    end
end
stimphaserel1 = stimphase1 + repmat(relphasemn,[1 1 opt.npre]);
[stimphaserel1,~,stimphaserel1std] = angmean(2*pi*flatten(stimphaserel1,2:3), 2);
stimphaserel1 = stimphaserel1/(2*pi);
stimphaserel1std = stimphaserel1std/(2*pi);

%and get the burst numbers relative to the pulses
[burstnum,burstpulse] = get_burst_pulse_num(data.burston,data.burstoff, ...
    pulset,pulsedur2);

if isfield(data,'stimcyclet');
    data = rmfield(data,{'stimcyclet','burstspercycle','burstcycle','burstcyclet'});
end

data.burstnum = burstnum;
data.burstpulse = burstpulse;

data.stimphase = stimphaserel1;
data.stimphasestd = stimphaserel1std;
data.pulset = pulset;
data.burstfreq = 1./burstper;
data.burstdur = data.burstoff - data.burston;

if opt.showdiagnostics    
    spiket1 = eventtrigger(data.pulset, data.spiket', [-4 5],'timerange');
    [burstt1,burstoff1,burstnum1] = eventtrigger(data.pulset, data.burstt', [-4 5], ...
        data.burstoff',burstnum', 'timerange');
    
    [~,ord] = sortrows([data.Direction data.stimphase]);

    % t0 = burstoff1(:,:,1);
    % t0 = t0(burstnum1(:,:,1) == -1);
    if data.amp > 0
        basefreq = repmat(data.stimfreq, size(spiket1));
    else
        basefreq = nanmean(flatten(1./burstperpre,2:3),2);
        basefreq = repmat(shiftdim(basefreq,-1), [size(spiket1,1) 1 nchan]);
    end
    t0 = data.pulset - stimphaserel1./basefreq(1,:,1)';
    
    spiket1 = spiket1 - repmat(t0', [size(spiket1,1) 1 nchan]);
    spiket1 = spiket1 .* basefreq; % + repmat(stimphaserel1',[size(spiket1,1) 1 nchan]);
    spiket1 = spiket1(:,ord,:);
    
    burstt1 = burstt1 - repmat(t0', [size(burstt1,1) 1 nchan]);
    burstt1 = burstt1 .* basefreq(1:size(burstt1,1),:,:); % + ...
        %repmat(stimphaserel1',[size(burstt1,1) 1 nchan]);
    burstoff1 = burstoff1 - repmat(t0', [size(burstt1,1) 1 nchan]);
    burstoff1 = burstoff1 .* basefreq(1:size(burstt1,1),:,:); % + ...
        %repmat(stimphaserel1',[size(burstt1,1) 1 nchan]);
    
    burstt1 = burstt1(:,ord,:);
    burstnum1 = burstnum1(:,ord,:);
    rep = repmat(1:length(pulset),[size(burstt1,1) 1 nchan]);
    
    stimend1 = stimphaserel1 + 0.5*pulsedur2.*basefreq(1,:,1)';
    
    figureseries('Raster');
    col = 'bgm';
    
    gc = find(data.goodchan);
    for i = 1:length(gc)
        c = gc(i);
        subplot(length(gc),1,i);
        
        raster(spiket1(:,:,c),'k');
        
        addplot(stimphaserel1(ord),1:length(ord),'r*-');
        addplot(stimend1(ord),1:length(ord),'r-');
        vertplot(-3:4,'y--');
        
        nd = mode(nduring(:,c));
        addplot(bursttpre(ord,c,end) - t0(ord) + (nduring(ord,c)+1)./basefreq(1,ord,c)', ...
            1:length(ord),'g-');
        addplot(bursttpre(ord,c,end) - t0(ord) + (nduring(ord,c)+1)./basefreq(1,ord,c)' + ...
            dphasepost(ord,c,1)./basefreq(1,ord,c)', ...
            1:length(ord),'g-');
        
%         bt2 = burstt1(:,:,c);
%         bn2 = burstnum1(:,:,c);
%         rep2 = rep(:,:,c);
%         for j = 0:2
%             isburst = bn2 == j;
%             addplot(bt2(isburst),rep2(isburst),[col(j+1) 'o-']);
%         end
    end
    
    figureseries('PRC');
    subplot(2,1,1);
    isleft = data.Direction == -1;
    plot(stimphaserel1(isleft),dphasepost(isleft,:,1),'ko');
    addplot(stimphaserel1(isleft),dphasepost(isleft,:,2),'r*');

    subplot(2,1,2);
    plot(stimphaserel1(isleft),dphasepost(~isleft,:,1),'ko');
    addplot(stimphaserel1(isleft),dphasepost(~isleft,:,2),'r*');
end


