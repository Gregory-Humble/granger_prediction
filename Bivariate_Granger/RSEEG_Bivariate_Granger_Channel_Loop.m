%% Bivariate looped over all channel combinations

%clears everything to start fresh and clean!
clear all; close all; clc;

cd('E:\Bivariate_Granger_Loop_Trial\Data');
workDir = 'E:\Bivariate_Granger_Loop_Trial';
saveDir = 'E:\Bivariate_Granger_Loop_Trial\Granger_Causality_RSEEG';
%mkdir(saveDir);
eeglab;

%% change .set file to .mat via EEGlab to prevent EEG struct to EEG obj errors
% IDlist = {'301', '302', '303', '305', '306', '307', '308', '309', '310', '311', '312', '313', '314', '315', '316', '317', '318', '319', '320', '322', '324', '325', '326', '327', '328', '329', '330', '331', '333', '335', '336', '338', '341', '343', '345', '346', '347', '348', '349', '351'};  
IDlist = {'302'};
% ConditionList = {'BL', 'END'}
ConditionList = {'BL'};

% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%  for t = 1:numel(IDlist);
%     for w = 1:numel(ConditionList);
% EEG = pop_loadset('filename', [IDlist{t} '_' ConditionList{w} '_7ICArejected.set'],'filepath','E:\\Bivariate_Granger_Loop_Trial\\Data\\');
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% save([IDlist{t} '_' ConditionList{w} '_7ICArejected.mat'])
%     end
%  end

%% create list of folders
cd(saveDir);
% for q = 1:numel(IDlist);
%     thisID = IDlist{q};
%     newDir = [saveDir, filesep, thisID];
%     mkdir(newDir);
% end
%% Load .mat file to run GC analysis
for thisID = 1:numel(IDlist);
   for thisCondition = 1:numel(ConditionList);
source = ['E:\Bivariate_Granger_Loop_Trial\Data', filesep, IDlist{thisID}, '_' ConditionList{thisCondition}, '_7ICArejected.mat']
load(source, 'ALLCOM', 'ALLEEG', 'CURRENTSET', 'CURRENTSTUDY', 'EEG', 'globalvars', 'LASTCOM', 'PLUGINLIST', 'STUDY');
%% change directory to output/save folder
cd(saveDir);
%%
%creates empty Granger Causality GC struct under EEG
EEG.GC = [];

%remove SO1 channel
EEG.allchan(43) = [];

% define autoregression parameters
order = 10;

%list of chan names
channames = ({EEG.allchan.labels});

%creating EEG.GC struct timetables consisting of 42 timetables for all channel names
vec = zeros([1,1000]); % vector of zeroes to initialise the timetable
TT = timetable(vec', 'SampleRate', 500); %create timetable with transposed vec zeroes into column format
EEG.GC = struct('AF3', TT, 'AF4', TT, 'F7', TT, 'F5', TT, 'F3', TT, 'F1', TT, 'FZ', TT, 'F2', TT, 'F4', TT, 'F6', TT, 'F8', TT, 'FC5', TT, 'FC3', TT, 'FC1', TT, 'FCZ', TT, 'FC2', TT, 'FC4', TT, 'FC6', TT, 'C5', TT, 'C3', TT, 'C1', TT, 'CZ', TT, 'C2', TT, 'C4', TT, 'C6', TT, 'M1', TT, 'M2', TT, 'P7', TT, 'P5', TT, 'P3', TT, 'P1', TT, 'PZ', TT, 'P2', TT, 'P4', TT, 'P6', TT, 'P8', TT, 'PO3', TT, 'POZ', TT, 'PO4', TT, 'O1', TT, 'OZ', TT, 'O2', TT);

% define channels to compute granger synchrony between
for aa = 1:numel(EEG.allchan)
    for bb = 1:numel(EEG.allchan)
chan1name = EEG.allchan(aa).labels
chan2name = EEG.allchan(bb).labels
        if bb == aa
            continue
        end

% find the index of those channels
chan1 = find( strcmpi(chan1name,{EEG.chanlocs.labels}) );
chan2 = find( strcmpi(chan2name,{EEG.chanlocs.labels}) );

% get AR coefficients and error from each signal
[Ax,Ex] = armorf(EEG.data(chan1,:,1),1,EEG.pnts,order);
[Ay,Ey] = armorf(EEG.data(chan2,:,1),1,EEG.pnts,order);

%%% reconstruct the data using the autoregressive coefficients
% x = zeros(1,EEG.pnts);
% y = zeros(1,EEG.pnts);
% 
% x(1:order) = EEG.data(chan1,1:order,1);
% y(1:order) = EEG.data(chan2,1:order,1);
% 
% 
% for i = order+1:EEG.pnts
%     
%     % initialize
%     thispointX = 0;
%     thispointY = 0;
%     
%     for ai=1:order
%         thispointX = thispointX + EEG.data(chan1,i-ai,1)*Ax(ai);
%         thispointY = thispointY + EEG.data(chan2,i-ai,1)*Ay(ai);
%     end
%     x(i-1) = thispointX;
%     y(i-1) = thispointY;
% end 

%% plot figure
% figure(1), clf
% subplot(211)
% plot(EEG.times,EEG.data(chan1,:,1),'b', EEG.times,x,'r')
% legend({'Real data';'Reconstructed from ARmodel'})
% 
% subplot(212)
% plot(EEG.times,EEG.data(chan2,:,1),'b', EEG.times,y,'r')
% legend({'Real data';'Reconstructed from ARmodel'})

%% Granger prediction

% Bivariate autoregression and associated error term
% [Axy,E] = armorf(EEG.data([chan1 chan2],:,1),1,EEG.pnts,order);
% 
% 
% % time-domain causal estimate
% granger_chan2_to_chan1 = log(Ex/E(1,1)); % Ex is the variance of the univariate AR errors - E is the variance of the bivariate AR errors
% granger_chan1_to_chan2 = log(Ey/E(2,2));
% 
% disp([ 'Granger prediction from ' chan1name ' to ' chan2name ' is ' num2str(granger_chan1_to_chan2) ]);
% disp([ 'Granger prediction from ' chan2name ' to ' chan1name ' is ' num2str(granger_chan2_to_chan1) ]);
% 
% % compute granger prediction over time

% initialize
    
x2yT = zeros(1,EEG.pnts);
y2xT = zeros(1,EEG.pnts);


% GC parameters
iwin   = 300; % window size in ms
iorder = 20;  % in ms


% convert window/order to points
win   = round(iwin/(1000/EEG.srate));
order = round(iorder/(1000/EEG.srate));

for timei=1:EEG.pnts-win
    
    % data from all trials in this time window
    % Data should be normalized before computing Granger estimates - to
    % make sure both channels are in the same scale. Both should be
    % stationary signals
    tempdata = zscore(reshape(EEG.data([chan1 chan2],timei:timei+win-1,1),2,win),0,2);
    
    % fit AR models (model estimation from bsmart toolbox)
    [Ax,Ex] = armorf(tempdata(1,:),1,win,order);
    [Ay,Ey] = armorf(tempdata(2,:),1,win,order);
    [Axy,E] = armorf(tempdata     ,1,win,order);
    
    % time-domain causal estimate
    
    y2xT(timei) = log(Ex/E(1,1));
    x2yT(timei) = log(Ey/E(2,2));

end
 %creates field names of specific channel comparisons based on the current
 %loop iteration.
 % the following line gives us the AR causal estimate for ONE time window
 % (model order in pnts)
 field1 = [chan1name, '_to_', chan2name];
 field2 = [chan2name, '_to_', chan1name];
 value1 = x2yT;
 value2 = y2xT;


% Creates a timetable with column vectors of the granger time series for
% each channel-channel  pair
EEG.GC.(chan1name) = addvars(EEG.GC.(chan1name), value1', 'NewVariableNames', {field1});
EEG.GC.(chan1name) = addvars(EEG.GC.(chan1name), value2', 'NewVariableNames', {field2});
  
        end
%% save timetable output     
thisDir = [saveDir,filesep,IDlist{thisID}]
cd(thisDir)
writetimetable(EEG.GC.(chan1name), [chan1name, '_', IDlist{thisID}, '_', ConditionList{thisCondition}, 'RSEEG_bivariate_granger_time_series_values.csv'])
%need to save EEG.GC as well
% need to save in individual subject folders 
    end
end
end

%% Save whole EEG struct
%save(['EEG_343_bivariate_granger_time_series_values.mat'], 'EEG')

%% Save whole workspace 
%save(['EEG_343_bivariate_granger_time_series_values.mat'])

%% Load whole EEG struct
%load(['EEG_343_bivariate_granger_time_series_values.mat'], 'EEG')

%% Load whole workspace 
%load(['EEG_343_bivariate_granger_time_series_values.mat'])

%% Plot figures
% plot single granger graph
% figure(1), clf, hold on
%     
% subplot(2,1,1)
% plot(EEG.times,)
% hold on
% plot(EEG.times,[fig2plotchan2],'r')
% legend({[ 'GC: ' [fig2plotchan1] ];[ 'GC: ' [fig2plotchan2] ]})
% 
% subplot(2,1,2)
% plot(EEG.times,dat_output.AF3_to_F3)
% hold on
% plot(EEG.times,dat_output.F3_to_AF3,'r')
% legend({[ 'GC: AF3 -> F3' ];[ 'GC: F3 -> AF3' ]})
% 
% title([ 'Window length: ' num2str(iwin) ' ms, order: ' num2str(iorder) ' ms' ])
% xlabel('Time (ms)')
% ylabel('Granger prediction estimate')
% set(gca,'xlim',[0 1700])
