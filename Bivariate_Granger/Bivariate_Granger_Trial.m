%% Bivariate looped over all channel combinations
clear all; close all; clc;

cd('D:\Bivariate_Granger_Loop_Trial');
workDir='D:\Bivariate_Granger_Loop_Trial\';
source = ['D:\Bivariate_Granger_Loop_Trial', filesep, '301_BL_7ICArejected.mat'];
load(source);
% define autoregression parameters
order = 14;
% define channels to compute granger synchrony between
EEG.allchan(43) = [] %remove SO1 channel

dat_output = []

for a = 1:numel(EEG.allchan)
    for b = 1:numel(EEG.allchan)
chan1name = EEG.allchan(a).labels
chan2name = EEG.allchan(b).labels
        if a == b
            continue
        end

% find the index of those channels
chan1 = find( strcmpi(chan1name,{EEG.chanlocs.labels}) )
chan2 = find( strcmpi(chan2name,{EEG.chanlocs.labels}) )

% get AR coefficients and error from each signal
[Ax,Ex] = armorf(EEG.data(chan1,:,1),1,EEG.pnts,order);
[Ay,Ey] = armorf(EEG.data(chan2,:,1),1,EEG.pnts,order);

%%% reconstruct the data using the autoregressive coefficients
x = zeros(1,EEG.pnts);
y = zeros(1,EEG.pnts);

x(1:order) = EEG.data(chan1,1:order,1);
y(1:order) = EEG.data(chan2,1:order,1);


for i = order+1:EEG.pnts
    
    % initialize
    thispointX = 0;
    thispointY = 0;
    
    for ai=1:order
        thispointX = thispointX + EEG.data(chan1,i-ai,1)*Ax(ai);
        thispointY = thispointY + EEG.data(chan2,i-ai,1)*Ay(ai);
    end
    x(i-1) = thispointX;
    y(i-1) = thispointY;
end 


% figure(1), clf
% subplot(211)
% plot(EEG.times,EEG.data(chan1,:,1),'b', EEG.times,x,'r')
% legend({'Real data';'Reconstructed from ARmodel'})
% 
% subplot(212)
% plot(EEG.times,EEG.data(chan2,:,1),'b', EEG.times,y,'r')
% legend({'Real data';'Reconstructed from ARmodel'})

% Granger prediction


% Bivariate autoregression and associated error term
[Axy,E] = armorf(EEG.data([chan1 chan2],:,1),1,EEG.pnts,order);


% time-domain causal estimate
granger_chan2_to_chan1 = log(Ex/E(1,1));
granger_chan1_to_chan2 = log(Ey/E(2,2));

disp([ 'Granger prediction from ' chan1name ' to ' chan2name ' is ' num2str(granger_chan1_to_chan2) ]);
disp([ 'Granger prediction from ' chan2name ' to ' chan1name ' is ' num2str(granger_chan2_to_chan1) ]);

% compute granger prediction over time

% initialize
    
x2yT = zeros(1,EEG.pnts);
y2xT = zeros(1,EEG.pnts);


% GC parameters
iwin   = 300; % window size in ms
iorder = 15;  % in ms


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
    dat_output(a,:) = y2xT
    dat_output(a,1001,:) = x2yT
end

    end
end

% draw lines
figure(15), clf, hold on

plot(EEG.times,x2yT)
plot(EEG.times,y2xT,'r')
legend({[ 'GC: ' chan1name ' -> ' chan2name ];[ 'GC: ' chan2name ' -> ' chan1name ]})

title([ 'Window length: ' num2str(iwin) ' ms, order: ' num2str(iorder) ' ms' ])
xlabel('Time (ms)')
ylabel('Granger prediction estimate')
set(gca,'xlim',[0 1700])