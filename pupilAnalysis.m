% Pupil analysis
% written by YH
% 6/11/2019

function D = pupilAnalysis
clear; clear all;

% input
pp = input('Plot trial data? (y/n) :  ', 's');

load('02MH13yh_eye.mat')

for trl=1:STARTtime(end,1)
    
    if pp == 'y'
        % Figure setting
        if mod(trl,20)==1
            figure(...
                'InvertHardcopy', 'off',...
                'Color', [1 1 1],...
                'Position', [0 0 1600 1600]);
        end
    end
    
    % Time
    c.stime = STARTtime(trl,2);
    c.etime = ENDtime(trl,2);
    c.dtime = MSGtime(MSGtime(:,1)==trl & MSGtime(:,3)==2,2); % Disk presentation
    
    % Rawdata(trl,time,x,y,pupil size, 0)
    c.trl_RawData = RawData(RawData(:,2)>=c.stime & RawData(:,2)<=c.etime,:);
    c.nData = length(c.trl_RawData);
    c.start_end = [c.trl_RawData(1,2)-c.dtime c.trl_RawData(end,2)-c.dtime];
    
    % Calculate max dilation speeds
    for ii = 2:c.nData-1
        c.maxDilationSpeeds(ii)...
            = max(abs((c.trl_RawData(ii,5)-c.trl_RawData(ii-1,5))...
            /(c.trl_RawData(ii,2)-c.trl_RawData(ii-1,2))),...
            abs((c.trl_RawData(ii+1,5)-c.trl_RawData(ii,5))...
            /(c.trl_RawData(ii+1,2)-c.trl_RawData(ii,2))));
    end
    c.mds = [NaN c.maxDilationSpeeds]';
    c.tRawData = [c.trl_RawData c.mds];
    
    % Calculate threshold
    n=1; %consistent value
    c.thresh = median(c.maxDilationSpeeds)...
        + (n * median(abs(c.maxDilationSpeeds-median(c.maxDilationSpeeds))));
    
    % Count blinks
    % EBLINK(trl,0,stime,etime,dur)
    c.trl_EBLINK = EBLINK(EBLINK(:,1)==trl,:);
    c.nBlink = size(c.trl_EBLINK,1);
    
    %% Detecting noise with blinks
    if c.nBlink>0
        % Change invalid data to NaN
        for ii = 1:c.nBlink
            c.tRawData...
                (c.tRawData(:,2)>=c.trl_EBLINK(ii,3)-200 ...
                & c.tRawData(:,2)<=c.trl_EBLINK(ii,4)+200 ...
                & c.tRawData(:,end)>c.thresh,5)...
                = NaN;
        end
    end
    
    %     %% Detecting noise without blinks
    %     c.tRawData...
    %         ( c.tRawData(:,end)>c.thresh,5)...
    %         = NaN;
    
    % Up sampling rate
    idx = c.tRawData(1,2)-c.dtime:c.tRawData(end,2)-c.dtime;
    c.ppldat(:,1) = idx;
    
    for ii = 1:length(idx)
        l = find(idx(ii)==c.tRawData(:,2)-c.dtime);
        if isempty(l)
            c.ppldat(ii,2) = NaN;
        else
            c.ppldat(ii,2) = c.tRawData(l,5);
        end
        clear ID
    end
    
    % Fill missing data
    % fillmiss_linear: linear interpolation
    % fillmiss: spline interpolation
    c.ppldat(:,2) = fillmiss_linear(c.ppldat(:,2));
    
    % Plot trial data
    if pp == 'y'
        if mod(trl,20)==0
            subplot(4,5,20)
        else
            subplot(4,5,mod(trl,20))
        end
        plot(c.trl_RawData(:,2)- c.dtime,c.trl_RawData(:,5),'r--'); hold on;
        %     plot(c.tRawData(:,2)- c.dtime,c.tRawData(:,5),'k');
        plot(c.ppldat(:,1),c.ppldat(:,2),'k')
        title(['Trial' num2str(trl)]);
    end
    
    % Normalize data
    c.dilation = (c.ppldat(c.ppldat(:,1)==0,2)-c.ppldat(1,2))/c.ppldat(1,2);
    c.ppldat(:,2) = (c.ppldat(:,2)-c.ppldat(1,2))/c.ppldat(1,2);
    
    D.trial(trl)=c;
    
    clear c idx
    
end

%% Calculate the mean of each condition

% Create summary data with -3000ms to 1000ms range
idx = -3000:1000;
allData = NaN(length(D.trial),length(idx));

for trl = 1:length(D.trial)
    allData(trl,(D.trial(trl).start_end(1)-idx(1)+1:D.trial(trl).start_end(1)-idx(1)+length(D.trial(trl).ppldat(:,2))))...
        = D.trial(trl).ppldat(:,2);
end

% Load behavior
load('02_MH_1_3.mat')
nn = 0;
for block = 1:length([Data.condition])
    for trl = 1:length(Data(block).TR)
        nn = nn + 1;
        tIdx(nn) = Data(block).TR(trl).train;
        cIdx(nn) = Data(block).condition;
        eIdx(nn) = Data(block).TR(trl).eyeError;
    end
end

% Mean
D.mean.single_c...
    = nanmean(allData(cIdx==1 & eIdx==0,:),1);
D.mean.single_p.trained...
    = nanmean(allData(tIdx==1 & cIdx==2 & eIdx==0,:),1);
D.mean.single_p.untrained...
    = nanmean(allData(tIdx==0 & cIdx==2 & eIdx==0,:),1);
D.mean.dual.trained...
    = nanmean(allData(tIdx==1 & cIdx==3 & eIdx==0,:),1);
D.mean.dual.untrained...
    = nanmean(allData(tIdx==0 & cIdx==3 & eIdx==0,:),1);


nn=0;
for cond = 1:3
    if cond == 1
        nn=nn+1;
        dilData=[D.trial(cIdx==cond & eIdx==0).dilation];
        [nOut,dd] = func_outlier(dilData(~isnan(dilData)));
        D.dilation(nn) = mean(dd);
        clear dilData nOut dd
    else
        for train = [1 0]
            nn=nn+1;
            dilData=[D.trial(tIdx==train & cIdx==cond & eIdx==0).dilation];
            [nOut,dd] = func_outlier(dilData(~isnan(dilData)));
            D.dilation(nn) = mean(dd);
            clear dilData nOut dd
        end
    end
end


% Plot with sdat to edat range
sdat = -3000;
edat = 1000;

figure(...
    'InvertHardcopy', 'off',...
    'Color', [1 1 1],...
    'Position', [0 0 1600 300]);

subplot(1,5,1);
plot(sdat:edat,D.mean.single_c(find(idx==sdat):find(idx==edat)));
title('Single central')
subplot(1,5,2);
plot(sdat:edat,D.mean.single_p.trained(find(idx==sdat):find(idx==edat)));
title('Single peripheral (trained)')
subplot(1,5,3);
plot(sdat:edat,D.mean.single_p.untrained(find(idx==sdat):find(idx==edat)));
title('Single peripheral (untrained)')
subplot(1,5,4);
plot(sdat:edat,D.mean.dual.trained(find(idx==sdat):find(idx==edat)));
title('Dual (trained)')
subplot(1,5,5);
plot(sdat:edat,D.mean.dual.untrained(find(idx==sdat):find(idx==edat)));
title('Dual (untrained)')

for ii=1:5
    subplot(1,5,ii)
    hold on;
    ylim([-.1 .4])
    plot([0 0],[-0.1 0.4],'--r');
end

return

function newdata = fillmiss(olddata)

% grab integer indices of v to input to spline, and find where NaNs are
v = olddata;
x = 1:length(v);
m = isnan(v);
% x(~m) contains the indices of the non-NaN values
% v(~m) contains the non-NaN values
% x(m) contains the indices of the NaN values, and thus the points at
% which we would like to query the spline interpolator
s = spline(x(~m),v(~m),x(m));

% replace NaN values with interpolated values; plot to see results
v(m) = s;
newdata = v;
return

function newdata = fillmiss_linear(olddata)

v = olddata;
x = 1:length(v);
m = isnan(v);

s =interp1(x(~m),v(~m),x,'linear');
newdata = s;
return

function [nOut,newdata] = func_outlier(olddata)
nOut=0; nData = length(olddata);
data = olddata;

% %% recursion
% 	while 1
% 		maxD = max(data);
% 		idx = find(data<maxD);
%
% 		if isempty(idx)
% 			break
% 		end
%
% 		mean_except_max = mean(data(idx));
% 		sd_except_max = std(data(idx));
%
% 		if maxD > mean_except_max + (sd_except_max*4);
% 			data = data(idx);
% 			nOut = nOut + 1;
% 		else
% 			break;
% 		end
% 	end
%
% % 	if nOut~=0
% % 		fprintf('excluded %d of %d\n',nOut,nData);
% % 	end
%
% 	newdata = data;

%% mean+nSD
idx=find(data<mean(data)+std(data)*3 & data>mean(data)-std(data)*3);
newdata = data(idx);
nOut = nData-length(newdata);

% %% cut under 1500ms
% idx=find(data<1.5);
% newdata = data(idx);
% nOut = nData-length(newdata);

return


