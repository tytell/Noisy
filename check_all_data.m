function check_all_data

basepath = '/Volumes/Data/CPG perturbations/'; %'E:\CPG perturbations';

outfile = 'perturbationsdata002.csv';

%files = {...

%      '2014-09-17/L29-rost-042.h5', ...
%     '2014-12-18/L33-caud-037.h5', ...
%     '2014-12-18/L33-caud-031.h5', ...      % not entrained
%     '2014-12-18/L33-caud-030.h5', ...
% };
 %'2014-07-11/L24-rost-052.h5', ...             % good file
%     '2014-07-11/L24-rost-059.h5', ...         % good file

inputyn('','clearsaved');

files = {...
%         '2016-04-15/L72-rost-014.tdms', ...
%         '2016-04-15/L72-rost-015.tdms', ...
%         '2016-04-15/L72-rost-019.tdms', ...
%         '2016-04-15/L72-rost-023.tdms',};
%      '2014-06-20/L22-rost-039.h5'    % L22
%      '2014-06-20/L22-rost-050.h5'
%      '2014-06-20/L22-rost-051.h5'
%      '2014-06-20/L22-rost-052.h5'
...%     '2014-07-11/L24-rost-052.h5'    % L24
...%     '2014-07-11/L24-rost-053.h5'
...%     '2014-07-11/L24-rost-054.h5'
...%     '2014-07-11/L24-rost-059.h5'
...%     '2014-07-11/L24-rost-059.h5'
...%     '2014-09-17/L29-rost-030.h5'    % L29
...%     '2014-09-17/L29-rost-031.h5'
...%     '2014-09-17/L29-rost-032.h5'
...%     '2014-09-17/L29-rost-042.h5'
...%BAD FILE    '2014-09-17/L29-rost-044.h5' 
...%     '2014-09-17/L29-rost-045.h5'
...%     '2014-09-17/L29-rost-046.h5'
...%     '2014-09-17/L29-rost-047.h5'
...%     '2014-09-17/L29-rost-048.h5'
...%%BAD FILE     '2014-10-17/L31-caud-018.h5'    % L31
...%BAD FILE     '2014-10-17/L31-caud-019.h5'
...%BAD FILE     '2014-10-17/L31-caud-020.h5'
...%     '2014-10-17/L31-caud-022.h5'
...%     '2014-10-17/L31-caud-023.h5'
...%     '2014-10-17/L31-caud-024.h5'
...%     '2014-10-17/L31-caud-025.h5'
...%     '2014-11-13/L32-caud-022.h5'    % L32-caud
...%%BAD FILE      '2014-11-13/L32-caud-023.h5'
%    '2014-11-13/L32-caud-027.h5'
%%BAD FILE      '2014-11-13/L32-caud-026.h5'
%    '2014-11-13/L32-caud-028.h5'
%    '2014-11-13/L32-caud-029.h5'
%    '2014-11-14/L32-rost-080.h5'    % L32-rost
%    '2014-11-14/L32-rost-081.h5'
%    '2014-11-14/L32-rost-083.h5'
%    '2014-11-14/L32-rost-084.h5'
%   '2014-11-14/L32-rost-085.h5'
 %    '2014-12-18/L33-caud-026.h5'    % L33-caud
%BAD FILE     '2014-12-18/L33-caud-027.h5'
%     '2014-12-18/L33-caud-028.h5'
%     '2014-12-18/L33-caud-029.h5'
%     '2014-12-18/L33-caud-030.h5'
%     '2014-12-18/L33-caud-031.h5'
%     '2014-12-18/L33-caud-037.h5'
%     '2014-12-18/L33-caud-038.h5'
%     '2014-12-18/L33-caud-039.h5'
%     '2014-12-19/L33-rost-063.h5'    % L33-rost
%     '2014-12-19/L33-rost-071.h5'
%      '2014-12-19/L33-rost-059.h5'
%     '2014-12-19/L33-rost-073.h5'
   %  '2015-02-11/L35-rost-021.h5'    % L35
   %  '2015-02-11/L35-rost-022.h5'
   %  '2015-02-11/L35-rost-023.h5'
   %  '2015-02-11/L35-rost-024.h5'
  %   '2015-02-11/L35-rost-025.h5'
%AD FILE     '2015-02-11/L35-rost-026.h5'
 %    '2015-02-11/L35-rost-027.h5'
 %    '2015-02-11/L35-rost-028.h5'
 %    '2015-02-11/L35-rost-033.h5'
 %    '2015-03-25/L40-rost-046.tdms'    % L40
 %    '2015-03-25/L40-rost-048.tdms'
 %   '2015-03-25/L40-rost-050.tdms'
 %    '2015-03-25/L40-rost-056.tdms'
 %   '2015-03-25/L40-rost-058.tdms'
 %    '2015-03-25/L40-rost-060.tdms'
 %    '2015-03-25/L40-rost-064.tdms'  %(0 amp)
 %    '2015-03-25/L40-rost-065.tdms'  %(0 amp)
 %    '2015-06-26/L48-caud-012.tdms'    % L48
 %    '2015-06-26/L48-caud-013.tdms'
 %    '2015-06-26/L48-caud-014.tdms'
 %    '2015-06-26/L48-caud-015.tdms'
  %   '2015-06-26/L48-caud-016.tdms'
  %   '2015-06-26/L48-caud-017.tdms'
  %   '2015-06-26/L48-caud-018.tdms'
  %   '2015-06-26/L48-caud-019.tdms'
 %    '2015-06-26/L48-caud-025.tdms'
%   %  '2015-06-26/L48-caud-026.tdms'
 %    '2015-06-26/L48-caud-027.tdms'
%     '2015-06-26/L48-caud-028.tdms'
  %   '2015-06-26/L48-caud-029.tdms'
 %    '2015-06-26/L48-caud-031.tdms'
%     '2015-06-26/L48-caud-032.tdms'
  %   '2015-06-26/L48-caud-033.tdms'
 %    '2015-06-26/L48-caud-038.tdms'
%     '2015-06-26/L48-caud-039.tdms'
 %    '2015-06-26/L48-caud-040.tdms'
%     '2015-06-26/L48-caud-041.tdms'
%     '2015-06-26/L48-caud-042.tdms'
%     '2015-06-26/L48-caud-043.tdms'
%    '2015-06-26/L48-caud-044.tdms'
%     '2015-06-26/L48-caud-045.tdms'
%BAD FILE     '2015-06-26/L48-caud-046.tdms'
%     '2015-06-26/L48-caud-048.tdms'
%     '2015-06-26/L48-caud-049.tdms'
%     '2015-06-26/L48-caud-050.tdms'
%BAD FILE     '2015-06-26/L48-caud-051.tdms'
%     '2015-06-26/L48-caud-052.tdms'
%     '2015-10-30/L59-rost-043.tdms'    % L59
%     '2015-10-30/L59-rost-044.tdms'
%     '2015-10-30/L59-rost-023.tdms'
%     '2015-10-30/L59-rost-024.tdms'
%     '2015-10-30/L59-rost-031.tdms'
%     '2015-10-30/L59-rost-032.tdms'
%     '2015-10-30/L59-rost-035.tdms'
%     '2015-10-30/L59-rost-036.tdms'
%     '2015-10-30/L59-rost-037.tdms'
%     '2015-10-30/L59-rost-038.tdms'
%     '2015-10-30/L59-rost-039.tdms'
%     '2015-10-30/L59-rost-040.tdms'
%     '2015-10-30/L59-rost-041.tdms'
%     '2015-10-30/L59-rost-042.tdms'
    '2015-11-13/L60-rost-018.tdms'    % L60
    '2015-11-13/L60-rost-014.tdms'
    '2015-11-13/L60-rost-015.tdms'
    '2015-11-13/L60-rost-018.tdms'
    '2015-11-13/L60-rost-019.tdms'
    '2015-11-13/L60-rost-022.tdms'
    '2015-11-13/L60-rost-023.tdms'
    '2015-11-13/L60-rost-024.tdms'
    '2015-11-13/L60-rost-025.tdms'
    '2015-11-13/L60-rost-027.tdms'
    '2015-11-13/L60-rost-029.tdms'
    '2015-11-13/L60-rost-030.tdms'
    '2015-11-13/L60-rost-031.tdms'
    '2015-11-13/L60-rost-032.tdms'
    '2015-11-13/L60-rost-033.tdms'
    '2015-11-13/L60-rost-034.tdms'
    };

stimuluslocations.L22_rost = 18; 
stimuluslocations.L24_rost = 24;        
stimuluslocations.L29_rost = 22;
stimuluslocations.L31_caud = 46;        
stimuluslocations.L32_caud = 47;        
stimuluslocations.L32_rost = 20;
stimuluslocations.L33_caud = 52;
stimuluslocations.L33_rost = 30;
stimuluslocations.L35_rost = 19;
stimuluslocations.L40_rost = 20;
stimuluslocations.L36_caud = 44;
stimuluslocations.L48_caud = 42; 
stimuluslocations.L59_rost = 22;
stimuluslocations.L60_rost = 21;
stimuluslocations.L72_rost = 18;


override_good_chan.L60 = [false true true true];
%override_good_chan.L60_rost_015 = [false true false true];
%override_good_chan.L60_caud = [false true false true];
override_good_chan.L40 = [false true true];


%from Daniela's processing code:
params.L24_rost_052.threshold = [-0.450816 -0.0402448 -0.524267 -0.277275;0.540617 0.0558435 0.339285 0.238219];
params.L24_rost_052.interburstdur = [0.12 0.3 0.13 0.12];
params.L24_rost_052.minspikes = [2 2 2 2];
params.L24_rost_052.goodchan = [1 0 1 1];
params.L24_rost_053.threshold = [-0.396072 -0.0665947 -0.396863 -0.240827;0.396072 0.0665947 0.396863 0.240827];
params.L24_rost_053.interburstdur = [0.18 0.3 0.15 0.15];
params.L24_rost_053.minspikes = [3 2 3 3];
params.L24_rost_053.goodchan = [1 0 1 1];
params.L24_rost_054.threshold = [-0.396072 -0.0665947 -0.552317 -0.260564;0.396072 0.0665947 0.396863 0.240827];
params.L24_rost_054.interburstdur = [0.18 0.3 0.15 0.15];
params.L24_rost_054.minspikes = [3 2 3 3];
params.L24_rost_054.goodchan = [1 0 1 1];
params.L24_rost_059.threshold = [-0.396072 -0.0665947 -0.706312 -0.259629;0.403157 0.0665947 0.291297 0.306765];
params.L24_rost_059.interburstdur = [0.12 0.3 0.15 0.12];
params.L24_rost_059.minspikes = [3 2 2 3];
params.L24_rost_059.goodchan = [1 0 1 1];
params.L29_rost_030.threshold = [-0.300613 -0.315143 -1.05835 -0.211557 -0.813337;0.211275 0.190755 0.829924 0.126043 0.813337];
params.L29_rost_030.interburstdur = [0.15 0.15 0.3 0.15 0.3];
params.L29_rost_030.minspikes = [2 2 2 2 2];
params.L29_rost_030.goodchan = [1 1 0 1 0];
params.L29_rost_031.threshold = [-0.37587 -0.258293 -0.130567 -0.216227 -0.295071;0.366466 0.204933 0.130567 0.137501 0.295071];
params.L29_rost_031.interburstdur = [0.2 0.15 0.3 0.12 0.3];
params.L29_rost_031.minspikes = [2 2 2 2 2];
params.L29_rost_031.goodchan = [1 1 0 1 0];
params.L29_rost_032.threshold = [-0.37587 -0.258293 -0.130567 -0.216227 -0.295071;0.366466 0.204933 0.130567 0.137501 0.295071];
params.L29_rost_032.interburstdur = [0.2 0.15 0.3 0.12 0.3];
params.L29_rost_032.minspikes = [2 2 2 2 2];
params.L29_rost_032.goodchan = [1 1 0 1 0];
params.L29_rost_042.threshold = [-0.300613 -0.315143 -1.05835 -0.211557 -0.813337;0.211275 0.190755 0.829924 0.126043 0.813337];
params.L29_rost_042.interburstdur = [0.15 0.15 0.3 0.15 0.3];
params.L29_rost_042.minspikes = [2 2 2 2 2];
params.L29_rost_042.goodchan = [1 1 0 1 0];
params.L29_rost_045.threshold = [-0.324747 -0.287135 -1.10633 -0.217001 -0.290525;0.346075 0.211035 0.996983 0.151408 0.290525];
params.L29_rost_045.interburstdur = [0.15 0.2 0.3 0.15 0.3];
params.L29_rost_045.minspikes = [2 2 2 2 2];
params.L29_rost_045.goodchan = [1 1 0 1 0];
params.L29_rost_046.threshold = [-0.300613 -0.315143 -1.05835 -0.211557 -0.813337;0.211275 0.190755 0.829924 0.126043 0.813337];
params.L29_rost_046.interburstdur = [0.15 0.15 0.3 0.15 0.3];
params.L29_rost_046.minspikes = [2 2 2 2 2];
params.L29_rost_047.threshold = [-0.278783 -0.257107 -0.370001 -0.122952 -0.183135;0.275804 0.2077 0.633096 0.168052 0.293775];
params.L29_rost_047.interburstdur = [0.15 0.2 0.3 0.2 0.3];
params.L29_rost_047.minspikes = [2 2 2 2 2];
params.L29_rost_047.goodchan = [1 1 0 1 0];
params.L32_caud_022.threshold = [-0.224697 -0.235924 -0.186034;0.120675 0.156165 0.222845];
params.L32_caud_022.interburstdur = [0.15 0.2 0.15];
params.L32_caud_022.minspikes = [2 2 2];
params.L32_caud_022.goodchan = [1 1 1];
params.L32_caud_027.threshold = [-0.241111 -0.100676 -0.127067;0.138337 0.111071 0.135523];
params.L32_caud_027.interburstdur = [0.3 0.3 0.15];
params.L32_caud_027.minspikes = [1 2 2];
params.L32_caud_027.goodchan = [1 1 1];
params.L32_caud_028.threshold = [-0.241111 -0.100676 -0.127067;0.138337 0.111071 0.135523];
params.L32_caud_028.interburstdur = [0.3 0.3 0.15];
params.L32_caud_028.minspikes = [1 2 2];
params.L32_caud_028.goodchan = [1 1 1];
params.L32_caud_023.threshold = [-0.241111 -0.125672 -0.195292;0.138337 0.138943 0.209799];
params.L32_caud_023.interburstdur = [0.3 0.3 0.3];
params.L32_caud_023.minspikes = [1 2 2];
params.L32_caud_023.goodchan = [1 1 1];
params.L32_caud_026.threshold = [-0.241111 -0.125672 -0.195292;0.138337 0.138943 0.209799];
params.L32_caud_026.interburstdur = [0.3 0.3 0.3];
params.L32_caud_026.minspikes = [1 2 2];
params.L32_caud_026.goodchan = [1 1 1];
params.L32_rost_081.threshold = [-0.151801 -0.274477 -0.193808;0.126155 0.231663 0.17973];
params.L32_rost_081.interburstdur = [0.2 0.2 0.2];
params.L32_rost_081.minspikes = [2 2 2];
params.L32_rost_081.goodchan = [1 1 1];
params.L32_rost_080.threshold = [-0.122331 -0.263708 -0.178244;0.123511 0.188587 0.164166];
params.L32_rost_080.interburstdur = [0.15 0.2 0.15];
params.L32_rost_080.minspikes = [2 2 2];
params.L32_rost_080.goodchan = [1 1 1];
params.L33_caud_026.threshold = [-0.218958 -0.315235 -0.341507;0.209078 0.315235 0.797589];
params.L33_caud_026.interburstdur = [0.2 0.15 0.15];
params.L33_caud_026.minspikes = [2 2 2];
params.L33_caud_026.goodchan = [1 1 1];
params.L33_caud_027.threshold = [-0.251432 -0.864881 -0.408651;0.224219 0.455656 0.564664];
params.L33_caud_027.interburstdur = [0.2 0.1 0.1];
params.L33_caud_027.minspikes = [2 2 2];
params.L33_caud_027.goodchan = [1 1 1];
params.L33_caud_028.threshold = [-0.218958 -0.59826 -0.341507;0.209078 0.401542 0.797589];
params.L33_caud_028.interburstdur = [0.15 0.15 0.15];
params.L33_caud_028.minspikes = [2 2 2];
params.L33_caud_028.goodchan = [1 1 1];
params.L33_caud_029.threshold = [-0.241813 -0.453836 -0.304586;0.265527 0.480014 0.285299];
params.L33_caud_029.interburstdur = [0.2 0.2 0.15];
params.L33_caud_029.minspikes = [2 2 2];
params.L33_caud_029.goodchan = [1 1 1];
params.L33_caud_030.threshold = [-0.241813 -0.453836 -0.304586;0.265527 0.480014 0.285299];
params.L33_caud_030.interburstdur = [0.15 0.2 0.12];
params.L33_caud_030.minspikes = [2 2 2];
params.L33_caud_030.goodchan = [1 1 1];
params.L33_caud_031.threshold = [-0.228542 -0.709174 -0.329383;0.235165 0.563859 0.595834];
params.L33_caud_031.interburstdur = [0.2 0.12 0.15];
params.L33_caud_031.minspikes = [2 2 2];
params.L33_caud_031.goodchan = [1 1 1];
params.L33_caud_037.threshold = [-0.228542 -0.709174 -0.329383;0.235165 0.479551 0.430921];
params.L33_caud_037.interburstdur = [0.2 0.12 0.15];
params.L33_caud_037.minspikes = [2 2 2];
params.L33_caud_037.goodchan = [1 1 1];
params.L33_rost_063.threshold = [-0.102531 -0.465189 -0.4087;0.139795 0.345513 0.430921];
params.L33_rost_063.interburstdur = [0.2 0.3 0.15];
params.L33_rost_063.minspikes = [2 2 2];
params.L33_rost_063.goodchan = [1 1 1];
params.L33_rost_071.threshold = [-0.102531 -0.401856 -0.4087;0.111304 0.345513 0.371518];
params.L33_rost_071.interburstdur = [0.2 0.2 0.15];
params.L33_rost_071.minspikes = [2 2 2];
params.L33_rost_071.goodchan = [1 1 1];
params.L35_rost_021.threshold = [-0.211056 -0.133545 -0.275463;0.195827 0.157625 0.309486];
params.L35_rost_021.interburstdur = [0.2 0.2 0.2];
params.L35_rost_021.minspikes = [1 2 1];
params.L35_rost_021.goodchan = [1 1 1];
params.L35_rost_022.threshold = [-0.211056 -0.133545 -0.275463;0.195827 0.157625 0.309486];
params.L35_rost_022.interburstdur = [0.2 0.2 0.2];
params.L35_rost_022.minspikes = [1 2 1];
params.L35_rost_022.goodchan = [1 1 1];
params.L35_rost_023.threshold = [-0.211056 -0.133545 -0.275463;0.195827 0.157625 0.309486];
params.L35_rost_023.interburstdur = [0.2 0.2 0.2];
params.L35_rost_023.minspikes = [1 2 1];
params.L35_rost_023.goodchan = [1 1 1];
params.L35_rost_024.threshold = [-0.202344 -0.144195 -0.275463;0.193162 0.146057 0.309486];
params.L35_rost_024.interburstdur = [0.2 0.15 0.2];
params.L35_rost_024.minspikes = [2 2 1];
params.L35_rost_024.goodchan = [1 1 1];
params.L35_rost_025.threshold = [-0.202344 -0.144195 -0.275463;0.193162 0.146057 0.309486];
params.L35_rost_025.interburstdur = [0.2 0.15 0.2];
params.L35_rost_025.minspikes = [2 2 1];
params.L35_rost_025.goodchan = [1 1 1];
params.L35_rost_033.threshold = [-0.202344 -0.144195 -0.275463;0.193162 0.146057 0.309486];
params.L35_rost_033.interburstdur = [0.2 0.15 0.2];
params.L35_rost_033.minspikes = [2 2 1];
params.L35_rost_033.goodchan = [1 1 1];

%check that all the paths are correct
good = true(size(files));
for i = 1:length(files)
    if ~exist(fullfile(basepath,files{i}),'file')
        good(i) = false;
    end
end
if any(~good)
    fprintf('Files not found:\n');
    fprintf('  %s\n', files{~good});
else
    fprintf('All files found.\n');
end

interburstdur = [];
threshold = [];
minspikes = [];
goodchan = [];

indiv = cell(size(files));
part = cell(size(files));
filenumstr = cell(size(files));
filenum = zeros(size(files));
stimloc = zeros(size(files));
for f = 1:length(files)
    fprintf('*** %s\n', files{f});
    [pn,fn,ext] = fileparts(files{f});
    
    tok = regexp(fn,'(L\d+)-(\w+)-(\d+)','tokens','once');
    assert(length(tok) == 3);
    
    indiv{f} = tok{1};
    part{f} = tok{2};
    filenumstr{f} = tok{3};
    filenum(f) = str2double(tok{3});
    
    if ~isfield(stimuluslocations, [indiv{f} '_' part{f}])
        error('Please add stimulus location for %s', [indiv{f} '_' part{f}]);
    end
    stimloc(f) = stimuluslocations.([indiv{f} '_' part{f}]);
    
    matfile = [fn '.mat'];
    if exist(matfile,'file') && ...
            inputyn('Load existing mat file? ','default',true)
        load(matfile,'data');
    else
        switch ext
            case '.h5'
                data = loadEntrain2Data(fullfile(basepath,files{f}), ...
                    'stimuluslocation',stimloc(f));
            case '.tdms'
                data = loadEntrain4Data(fullfile(basepath,files{f}), ...
                    'stimuluslocation',stimloc(f));
        end
        
        nchan = size(data.sig,2);

        fnvar = fn;
        fnvar(fn == '-') = '_';
        if isfield(params,fnvar)
            interburstdur = params.(fnvar).interburstdur;
            threshold = params.(fnvar).threshold;
            minspikes = params.(fnvar).minspikes;
            
            if isfield(params.(fnvar),'goodchan')
                goodchan = params.(fnvar).goodchan;
            else
                goodchan = true(1,nchan);
            end
        else
            interburstdur = [];
            threshold = [];
            minspikes = [];
            goodchan = [];
        end
        
        pulsetimes = data.Time+data.BeforeDurSec+data.Phase/data.stimfreq;
        data = findbursts_gui(data, 'threshold',threshold, ...
            'interburstdur',interburstdur, 'minspikes',minspikes, ...
            'goodchan',goodchan, 'eventtimes',pulsetimes);
    end
    
    if ~isfield(data,'bursttpre') || inputyn('Process pulse data again?','default',false)
        data = process_pulse_data2(data, 'savediagnostics','diagnosticname',[fn '.pdf']);
        
        save(matfile,'data');
    end
    if sum(abs(data.dphasepost(:)) > 0.5) > 0.1 * sum(isfinite(data.dphasepost(:)))
        warning('Phase matching seems to be off');
    end
    
    goodchan = data.goodchan;
    if isfield(override_good_chan, [indiv{f} '_' part{f} '_' filenumstr{f}])
        goodchan = override_good_chan.([indiv{f} '_' part{f} '_' filenumstr{f}]);
    elseif isfield(override_good_chan, [indiv{f} '_' part{f}])
        goodchan = override_good_chan.([indiv{f} '_' part{f}]);
    elseif isfield(override_good_chan, indiv{f})
        goodchan = override_good_chan.(indiv{f});
    end
    if length(goodchan) ~= size(data.sig,2)
        warning('Override good channels should have length %d, but actually has length %d. Skipping', ...
            size(data.sig,2), length(goodchan));
        goodchan = data.goodchan;
    end
    nchan = sum(goodchan);
    
    nbursts = size(data.bursttpre,3) + size(data.bursttpost,3) + 1;
    npulses = size(data.bursttpre,1);

    t0 = data.bursttpre(:,data.stablechan,end);

    out.filename = repmat({fn},[npulses nchan nbursts]);
    out.indiv = repmat(indiv(f),[npulses nchan nbursts]);
    out.filenum = repmat(filenum(f), [npulses nchan nbursts]);
    out.channel = repmat(makerow(data.channelnames(goodchan)), [npulses 1 nbursts]);
    
    out.stimfreq = repmat(data.stimfreq, [npulses nchan nbursts]);
    out.stimamp = repmat(data.amp,[npulses nchan nbursts]);
    
    out.pulsephase = repmat(data.stimphase,[1 nchan nbursts]);
    out.pulsephasestd = repmat(data.stimphasestd,[1 nchan nbursts]);
    out.pulsedir = repmat(data.Direction,[1 nchan nbursts]);
    out.pulset = repmat(data.pulset - t0, [1 nchan nbursts]);
    out.origphase = repmat(data.Phase, [1 nchan nbursts]);                       
    
    out.burstnum = repmat(shiftdim(-size(data.bursttpre,3):size(data.bursttpost,3),-1), ...
        [npulses nchan 1]);
    
    out.burstt = cat(3,data.bursttpre,data.bursttduring,data.bursttpost);
    out.burstt = out.burstt(:,goodchan,:) - repmat(t0,[1 nchan nbursts]);
    out.burstfreq = 1./cat(3,data.burstperpre,data.burstperduring,data.burstperpost);
    out.burstfreq = out.burstfreq(:,goodchan,:);
    
    burstfreqpre = 1./nanmedian(flatten(data.burstperpre(:,goodchan,:),2:3),2);
    out.burstfreqpre = repmat(burstfreqpre,[1,nchan,nbursts]);
    
    out.dburstfreq = out.burstfreq - out.burstfreqpre;
    
    out.dphase = cat(3,data.dphasepre(:,goodchan,:),...
        NaN(npulses,nchan,2),...
        data.dphasepost(:,goodchan,:));
    out.donphase = cat(3,data.donphasepre(:,goodchan,:),...
        NaN(npulses,nchan,2),...
        data.donphasepost(:,goodchan,:));
    out.doffphase = cat(3,data.doffphasepre(:,goodchan,:),...
        NaN(npulses,nchan,2),...
        data.doffphasepost(:,goodchan,:));
    
    out.relburstdur = cat(3,data.burstdurfracpre,data.burstdurfracduring,...
        data.burstdurfracpost);
    out.relburstdur = out.relburstdur(:,goodchan,:);
    
    fields = fieldnames(out);
    for i = 1:length(fields)
        fieldname = fields{i};
        out.(fieldname) = permute(out.(fieldname), [2 3 1]);
    end
    
    save_struct_as_table(outfile,out, 'append',f > 1);
end
            
    