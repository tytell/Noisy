function [pers,freq,phase] = est_burst_period(burstt, t, varargin)

opt.nsmoothper = 3;
opt.discard = logical([]);

opt = parsevarargin(opt, varargin, 2);

nchan = size(burstt,2);
nburst = size(burstt,1);

per0 = NaN(size(burstt));
per0(1:end-1,:) = diff(burstt);
good = ~isnan(burstt);

isdouble = false(size(burstt));
isskip = false(size(burstt));

[t1,ord] = sort(burstt(good));
per1 = per0(good);
per1 = per1(ord);

if ~isempty(opt.discard)
    per1(opt.discard) = NaN;
end

pers1 = smooth(t1,per1, opt.nsmoothper*nchan, 'rlowess');

pers2 = NaN(size(per1));
pers2(ord) = pers1;

pers = NaN(size(per0));
pers(good) = pers2;

for c = 1:nchan
    tc = burstt(:,c);
    tc(~good(:,c)) = NaN;

    i = 1;
    while i <= nburst-1
        if ~good(i,c) || isnan(tc(i))
            i = i+1;
            continue;
        end

        %first collect bursts within 1.5 periods, skipping NaNs
        off = 0;
        while (i+off+1 <= nburst) && ~(tc(i+off+1) - tc(i) > 1.5*pers(i,c))
            off = off+1;
        end

        if off == 0
            %look for a skipped burst
            nskip = round((tc(i+1) - tc(i)) / pers(i,c));
            per0(i,c) = per0(i,c) / nskip;
            isskip(i,c) = true;
            i = i+1;
        elseif off > 1
            nper = (burstt(i+(1:off),c) - burstt(i,c)) / pers(i,c);
            isdouble1 = true(off,1);
            [~,ind] = min(abs(nper - 1));
            isdouble1(ind) = false;

            good(i+(1:off),c) = ~isdouble1;
            isdouble(i+(1:off),c) = isdouble1;

            i = i+ind;
        else
            i = i+1;
        end
    end
end

[t1,ord] = sort(burstt(good));
per1 = per0(good);
per1 = per1(ord);

if ~isempty(opt.discard)
    per1(opt.discard) = NaN;
end

pers1 = smooth(t1,per1, opt.nsmoothper*nchan, 'rlowess');

pers2 = NaN(size(per1));
pers2(ord) = pers1;

pers = NaN(size(per0));
pers(good) = pers2;

span = (t >= min(t1)) & (t <= max(t1));
good = isfinite(t1) & isfinite(pers1);
freq = NaN(size(t));
freq(span) = interp1(t1(good),1./pers1(good), t(span), 'cubic');
phase = NaN(size(t));
phase(span) = cumtrapz(t(span),freq(span));







