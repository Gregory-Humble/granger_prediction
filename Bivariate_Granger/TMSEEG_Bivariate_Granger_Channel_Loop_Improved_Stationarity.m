%% Bivariate Granger for TMS-EEG data
%% Load data and subtract trial average for each channel and time point

% clears everything to start fresh and clean!
clear all; close all; clc;

cd('E:\Bivariate_Granger_Loop_Trial\Data_TMS_EEG\trial_average_removed');
workDir = 'E:\Bivariate_Granger_Loop_Trial';
saveDir = 'E:\Bivariate_Granger_Loop_Trial\Granger_Causality_TMS_EEG\Model_Order_30';

eeglab;

%% change .set file to .mat via EEGlab to prevent EEG struct to EEG obj errors
% 305 (no rightpfc)
% 309 (no rightpfc)
% 316 (no leftpfc)
% 349 (no rightpfc)

% IDlist = {'301', '302', '305', '306', '307', '308', '309', '310', '311', '312', '314', '315', '316', '317', '318', '319', '320', '322', '324', '325', '326', '327', '328', '329', '330', '331', '333', '335', '336', '338', '341', '343', '345', '346', '347', '348', '349', '351'};  
IDlist = {'301', '302', '305', '306', '307', '308', '309', '310', '311', '312', '314', '315', '317', '318', '319', '320', '322', '324', '325', '326', '327', '328', '329', '330', '333', '335', '336', '338', '341', '343', '345', '346', '347', '348', '351'};  

% Regions = {'leftpfc', 'rightpfc'};
Regions = {'leftpfc'};

  for thisID = 1:numel(IDlist);
     for Reg = 1:numel(Regions);
EEG = pop_loadset('filename', [IDlist{thisID} '_TMSEEG_BL_' Regions{Reg} '_ds_ica1_filt_ica2_clean_reref.set'],'filepath','E:\Alz_Clinical_Trial\Alz_Data_Analysis_10JAN20\TMSEEG\TMSEEG_ALZ_clean_2\BL');
EEG = eeg_checkset(EEG);
save([IDlist{thisID} '_TMSEEG_BL_' Regions{Reg} '_ds_ica1_filt_ica2_clean_reref.mat'], 'EEG');
     end
  end
%% remove mean trial average (TEP) from each time point and channel
% cd('E:\Bivariate_Granger_Loop_Trial\Data_TMS_EEG');
%   for thisID = 1:numel(IDlist);
%     for Reg = 1:numel(Regions);
% cd('E:\Bivariate_Granger_Loop_Trial\Data_TMS_EEG');
% source = ['E:\Bivariate_Granger_Loop_Trial\Data_TMS_EEG', filesep, IDlist{thisID}, '_TMSEEG_BL_' Regions{Reg}, '_ds_ica1_filt_ica2_clean_reref.mat'];
% load(source, 'EEG')
%  for chani = 1:numel(EEG.nbchan);
%      for ti = 1:numel(EEG.pnts);
%      EEG.data(chani,ti,:) = EEG.data(chani,ti,:) - mean(EEG.data(chani,ti,:)); 
%      end
%  end
% cd('E:\Bivariate_Granger_Loop_Trial\Data_TMS_EEG\trial_average_removed');
% save([IDlist{thisID} '_TMSEEG_BL_' Regions{Reg} '_ds_ica1_filt_ica2_clean_reref_trialavgremoved.mat'], 'EEG');
%      end
%  end

%% create list of folders
cd(saveDir);
for foldList = 1:numel(IDlist);
     thisID = IDlist{foldList};
     newDir = [saveDir, filesep, thisID];
    mkdir(newDir);
end
cd('E:\Bivariate_Granger_Loop_Trial\Data_TMS_EEG\trial_average_removed');
%% Load .mat file to run GC analysis
for thisID = 1:numel(IDlist);
   for Reg = 1:numel(Regions);
source = ['E:\Bivariate_Granger_Loop_Trial\Data_TMS_EEG\trial_average_removed', filesep, IDlist{thisID}, '_TMSEEG_BL_' Regions{Reg}, '_ds_ica1_filt_ica2_clean_reref_trialavgremoved.mat']
load(source);

%% remove SO1, M1 and M2 channels
EEG.NoCh = {'FP1'; 'FPZ'; 'FP2'; 'FT7'; 'FT8'; 'T7'; 'T8'; 'TP7'; 'CP5'; 'CP3'; 'CP1'; 'CPZ'; 'CP2'; 'CP4'; 'CP6'; 'TP8'; 'PO7'; 'PO5'; 'PO6'; 'PO8'; 'CB1'; 'CB2'; 'E3'; 'HEOG'; 'M1'; 'M2'; 'SO1'}; 
EEG = pop_select(EEG,'nochannel',EEG.NoCh); 
EEG.allchan=EEG.chanlocs;

%% resize data to 0-1 seconds
EEG.data = EEG.data(:,1001:2001,:); % takes all channels, EEG.pnts between 1001-2001, and all trials
EEG.pnts = length(EEG.data); % updates the pnts variable
EEG.xmin = 0; % updates xmin
EEG.xmax = 1; % updates xmax
EEG.times = EEG.times(1001:2001); % updates the times vector
%% change directory to output/save folder
cd(saveDir);
%%
%creates empty Granger Causality GC struct under EEG
EEG.GC = [];

% define autoregression parameters
order = 30;

%list of chan names
channames = ({EEG.allchan.labels});

vec = zeros([1,1001]); % vector of zeroes to initialise the timetable
TT = timetable(vec', 'SampleRate', 1000); %create timetable with transposed vec zeroes into column format
%EEG.GC = struct('AF3', TT, 'AF4', TT, 'F7', TT, 'F5', TT, 'F3', TT, 'F1', TT, 'FZ', TT, 'F2', TT, 'F4', TT, 'F6', TT, 'F8', TT, 'FC5', TT, 'FC3', TT, 'FC1', TT, 'FCZ', TT, 'FC2', TT, 'FC4', TT, 'FC6', TT, 'C5', TT, 'C3', TT, 'C1', TT, 'CZ', TT, 'C2', TT, 'C4', TT, 'C6', TT, 'M1', TT, 'M2', TT, 'P7', TT, 'P5', TT, 'P3', TT, 'P1', TT, 'PZ', TT, 'P2', TT, 'P4', TT, 'P6', TT, 'P8', TT, 'PO3', TT, 'POZ', TT, 'PO4', TT, 'O1', TT, 'OZ', TT, 'O2', TT);
%for analysis of sepcifc electrodes - F3(, F4, P3, P4.
EEG.GC = struct('F3', TT, 'F4', TT, 'P3', TT, 'P4', TT);

% define channels to compute granger synchrony between - F3 (5), F4(9),
% P3(28), P4(32)
specificchans = [5, 9, 28, 32];
for aa = 1:numel(specificchans)
    for bb = 1:numel(specificchans)
currentchan1 = specificchans(aa);        
currentchan2 = specificchans(bb);  
chan1name = EEG.allchan(currentchan1).labels
chan2name = EEG.allchan(currentchan2).labels
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
iorder = 30;  % in ms


% convert window/order to points
win   = round(iwin/(1000/EEG.srate));
order = round(iorder/(1000/EEG.srate));
EEG.order = order % for saving model order in EEG struct if needed

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
cd(thisDir);
writetimetable(EEG.GC.(chan1name), [chan1name, '_', IDlist{thisID}, '_', Regions{Reg}, '_TMSEEG_bivariate_granger_time_series_values.csv'])
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
