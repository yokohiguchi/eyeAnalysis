% Pupil analysis
% written by YH
% 6/11/2019

function D = pupilAnalysis
clear; clear all;

% input
gType = input('Group type? (1:with training, 2:w/o training) :  ');
tType = input('Test type? (1:pre, 2:post, 3:re) :  ');
pp = input('Plot trial data? (y/n) :  ', 's');

if gType == 1
    ID = ['02';'03';'04';'05';'06';'08';'09';'11';'12';'14';'15';'16';'17';'21';'22';'24'];
    name = ['MH';'JY';'MN';'YM';'AU';'TM';'KK';'NW';'TW';'KT';'AS';'HK';'TI';'MK';'NK';'NF'];
    if tType == 1
        el_dir = 'data/eye_pre/';
        behav_dir = 'data/mat_pre/';
        fig_dir = 'figure/pupil_pre/';
    elseif tType == 2
        el_dir = 'data/eye_post/';
        behav_dir = 'data/mat_post/';
        fig_dir = 'figure/pupil_post/';
    elseif tType == 3
        el_dir = 'data/eye_re/';
        behav_dir = 'data/mat_re/';
        fig_dir = 'figure/pupil_re/';
    end
else
    ID = ['01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16'];
    name = ['RC';'MT';'KO';'EI';'MY';'MM';'SM';'RU';'YU';'NW';'RF';'RK';'KM';'MM';'TM';'DS'];
    if tType == 1
        el_dir = 'data/eye_pre_ctrl/';
        behav_dir = 'data/mat_pre_ctrl/';
        fig_dir = 'figure/pupil_pre_ctrl/';
    elseif tType == 2
        el_dir = 'data/eye_post_ctrl/';
        behav_dir = 'data/mat_post_ctrl/';
        fig_dir = 'figure/pupil_post_ctrl/';
    end
    
end

nSubject = length(ID);

for sub = 1:nSubject
    
    fprintf('.');
    
    for bb = 1:4
        
        % Load data
        if tType == 1
            if bb == 1
                sess = 1; r = 3;
            else
                sess = 2; r = bb-1;
            end
            el...
                = strcat(el_dir,ID(sub,:),name(sub,:),...
                num2str(sess),num2str(r),'yh_eye.mat');
            behav...
                = strcat(behav_dir,ID(sub,:),'_',name(sub,:),...
                '_',num2str(sess),'_',num2str(r),'.mat');
        else
            el...
                = strcat(el_dir,ID(sub,:),name(sub,:),...
                'T',num2str(bb),'yh_eye.mat');
            behav...
                = strcat(behav_dir,ID(sub,:),'_',name(sub,:),...
                '_T_',num2str(bb),'.mat');
        end
        
        load(el);
        load(behav);
        
        nn=0;empIdx = [];
        for trl=1:MSGtime(end,1)
            
            if pp == 'y'
                % Figure setting
                if mod(trl,20)==1
                    nn=nn+1;
                    h = figure(...
                        'InvertHardcopy', 'off',...
                        'Color', [1 1 1],...
                        'Position', [0 0 1600 1600]);
                end
            end
            
            % Time
            c.stime = STARTtime(trl,2);
            c.etime = ENDtime(trl,2);
            c.dtime = MSGtime(MSGtime(:,1)==trl & MSGtime(:,3)==2,2); % Disk presentation
            
            if isempty(c.dtime)
                c.dtime = MSGtime(MSGtime(:,1)==trl & MSGtime(:,3)==1,2)+8;
            end
            
            % Rawdata(trl,time,x,y,pupil size, 0)
            c.trl_RawData = RawData(RawData(:,2)>=c.stime & RawData(:,2)<=c.etime,:);
            c.nData = size(c.trl_RawData,1);
            
            if c.nData < 50
                c.start_end = [NaN NaN];
                c.tRawData = NaN;
                c.maxDilationSpeeds = NaN;
                c.mds = NaN;
                c.thresh = NaN;
                c.trl_EBLINK = [];
                c.nBlink = NaN;
                c.ppldat = NaN;
                c.dilation = NaN;
                D.sub(sub).data(bb).trial(trl)=c;
                empIdx = [empIdx trl];
                clear c idx
            else
                if gType == 2
                    c.trl_RawData(:,5) =(c.trl_RawData(:,5)/2).^2*pi/10000;
                end
                
                c.start_end = [c.trl_RawData(1,2)-c.dtime c.trl_RawData(end,2)-c.dtime];
                
                c.tRawData = c.trl_RawData;
                
                %                 % Detect noise
                %                 if ~isempty(c.tRawData(c.tRawData(1,5)-c.tRawData(:,5)>300,5))
                %                     c.tRawData(c.tRawData(1,5)-c.tRawData(:,5)>300,:)=[];
                %                 end
                
                % Calculate max dilation speeds
                for ii = 2:length(c.tRawData)-1
                    c.maxDilationSpeeds(ii)...
                        = max(abs((c.tRawData(ii,5)-c.tRawData(ii-1,5))...
                        /(c.tRawData(ii,2)-c.tRawData(ii-1,2))),...
                        abs((c.tRawData(ii+1,5)-c.tRawData(ii,5))...
                        /(c.tRawData(ii+1,2)-c.tRawData(ii,2))));
                end
                c.mds = [NaN c.maxDilationSpeeds]';
                c.tRawData = [c.tRawData c.mds];
                
                % Calculate threshold
                n=1; %consistent value
                c.thresh = median(c.maxDilationSpeeds)...
                    + (n * median(abs(c.maxDilationSpeeds-median(c.maxDilationSpeeds))));
                
                % Count blinks
                % EBLINK(trl,0,stime,etime,dur)
                c.trl_EBLINK = EBLINK(EBLINK(:,1)==trl,:);
                c.nBlink = size(c.trl_EBLINK,1);
                
                % Detect noise with blinks
                if c.nBlink>0
                    % Change invalid data to NaN
                    for ii = 1:c.nBlink
                        c.tRawData...
                            (c.tRawData(:,2)>=c.trl_EBLINK(ii,3) ...
                            & c.tRawData(:,2)<=c.trl_EBLINK(ii,4)+200,5)= NaN;
                        c.tRawData...
                            (c.tRawData(:,2)<=c.trl_EBLINK(ii,3) ...
                            & c.tRawData(:,2)>=c.trl_EBLINK(ii,3)-30,5)= NaN;
                    end
                end
                
                % Detect noise
                c.tRawData...
                    (c.tRawData(:,end)>c.thresh,5)...
                    = NaN;
                
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
                    plot(c.tRawData(:,2)- c.dtime,c.tRawData(:,5),'k');
                    plot(c.ppldat(:,1),c.ppldat(:,2),'k')
                    yl=get(gca,'ylim');
                    xl=get(gca,'xlim');
                    text(xl(1),yl(2),[num2str(c.nBlink) 'blinks']);
                    title(['Trial' num2str(trl)]);
                    drawnow;
                    
                    if mod(trl,20)==0 ||trl==MSGtime(end,1)
                        
                        if tType == 1
                            filename...
                                = strcat(fig_dir,ID(sub,:),name(sub,:),...
                                num2str(sess),num2str(r),'_',num2str(nn),'.fig');
                        else
                            filename...
                                = strcat(fig_dir,ID(sub,:),name(sub,:),...
                                'T',num2str(bb),'_',num2str(nn),'.fig');
                        end
                        
                        savefig(h,filename);
                        close(h);
                    end
                end
                
                % Normalize data
                c.dilation = (c.ppldat(c.ppldat(:,1)==0,2)-c.ppldat(1,2))/c.ppldat(1,2);
                c.ppldat(:,2) = (c.ppldat(:,2)-c.ppldat(1,2))/c.ppldat(1,2);
                
                D.sub(sub).data(bb).trial(trl)=c;
            end % if c.nData == 0
            clear c idx
        end % for trl
        
        % Create summary data with -3000ms to 1000ms range
        idx = -3000:1000; dd = D.sub(sub).data(bb);
        allData = NaN(length(dd.trial),length(idx));
        
        for trl = 1:length(dd.trial)
            if dd.trial(trl).start_end(1)<idx(1)
                continue;
            elseif ismember(trl,empIdx)
                continue;
            else
                allData(trl,(dd.trial(trl).start_end(1)-idx(1)+1:dd.trial(trl).start_end(1)-idx(1)+length(dd.trial(trl).ppldat(:,2))))...
                    = dd.trial(trl).ppldat(:,2);
            end
        end
        
        % Load behavior
        nn = 0;
        for block = 1:length([Data.condition])
            for trl = 1:length(Data(block).TR)
                nn = nn + 1;
                tIdx(nn) = Data(block).TR(trl).train;
                cIdx(nn) = Data(block).condition;
                eIdx(nn) = Data(block).TR(trl).eyeError;
            end
        end
        
        % Calculate the mean of each block
        D.sub(sub).single_c(bb,:)...
            = nanmean(allData(cIdx==1 & eIdx==0,:),1);
        D.sub(sub).single_p.trained(bb,:)...
            = nanmean(allData(tIdx==1 & cIdx==2 & eIdx==0,:),1);
        D.sub(sub).single_p.untrained(bb,:)...
            = nanmean(allData(tIdx==0 & cIdx==2 & eIdx==0,:),1);
        D.sub(sub).dual.trained(bb,:)...
            = nanmean(allData(tIdx==1 & cIdx==3 & eIdx==0,:),1);
        D.sub(sub).dual.untrained(bb,:)...
            = nanmean(allData(tIdx==0 & cIdx==3 & eIdx==0,:),1);
        
        
        nn=0;
        for cond = 1:3
            if cond == 1
                nn=nn+1;
                dilData=[dd.trial(cIdx==cond & eIdx==0).dilation];
                [nOut,newdata] = func_outlier(dilData(~isnan(dilData)));
                D.sub(sub).dilation(bb,nn) = mean(newdata);
                clear dilData nOut newdata
            else
                for train = [1 0]
                    nn=nn+1;
                    dilData=[dd.trial(tIdx==train & cIdx==cond & eIdx==0).dilation];
                    [nOut,newdata] = func_outlier(dilData(~isnan(dilData)));
                    D.sub(sub).dilation(bb,nn) = mean(newdata);
                    clear dilData nOut newdata
                end
            end
        end
        
        clear allData *Idx
    end % for bb
    
    % Calculate the mean of each condition
    D.mean.single_c(sub,:) = nanmean(D.sub(sub).single_c);
    D.mean.single_p.trained(sub,:) = nanmean(D.sub(sub).single_p.trained);
    D.mean.single_p.untrained(sub,:) = nanmean(D.sub(sub).single_p.untrained);
    D.mean.dual.trained(sub,:) = nanmean(D.sub(sub).dual.trained);
    D.mean.dual.untrained(sub,:) = nanmean(D.sub(sub).dual.untrained);
    D.mean.dilation(sub,:) = nanmean(D.sub(sub).dilation);
    
end % for sub

% Plot with sdat to edat range
sdat = -3000;
edat = 1000;

figure(...
    'InvertHardcopy', 'off',...
    'Color', [1 1 1],...
    'Position', [0 0 1600 600]);

subplot(2,5,1);
plot(sdat:edat,D.mean.single_c(:,find(idx==sdat):find(idx==edat))')
title('Single central')
subplot(2,5,2);
plot(sdat:edat,D.mean.single_p.trained(:,find(idx==sdat):find(idx==edat))');
title('Single peripheral (trained)')
subplot(2,5,3);
plot(sdat:edat,D.mean.single_p.untrained(:,find(idx==sdat):find(idx==edat))');
title('Single peripheral (untrained)')
subplot(2,5,4);
plot(sdat:edat,D.mean.dual.trained(:,find(idx==sdat):find(idx==edat))');
title('Dual (trained)')
subplot(2,5,5);
plot(sdat:edat,D.mean.dual.untrained(:,find(idx==sdat):find(idx==edat))');
title('Dual (untrained)')

subplot(2,5,6);
plot(sdat:edat,...
    nanmean(D.mean.single_c(:,find(idx==sdat):find(idx==edat)))',...
    'k','LineWidth',2);
title('Single central')
subplot(2,5,7);
plot(sdat:edat,...
    nanmean(D.mean.single_p.trained(:,find(idx==sdat):find(idx==edat)))',...
    'k','LineWidth',2);
title('Single peripheral (trained)')
subplot(2,5,8);
plot(sdat:edat,nanmean(D.mean.single_p.untrained(:,find(idx==sdat):find(idx==edat)))',...
    'k','LineWidth',2);
title('Single peripheral (untrained)')
subplot(2,5,9);
plot(sdat:edat,nanmean(D.mean.dual.trained(:,find(idx==sdat):find(idx==edat)))',...
    'k','LineWidth',2);
title('Dual (trained)')
subplot(2,5,10);
plot(sdat:edat,nanmean(D.mean.dual.untrained(:,find(idx==sdat):find(idx==edat)))',...
    'k','LineWidth',2);
title('Dual (untrained)')

for ii=1:10
    subplot(2,5,ii)
    hold on;
    
    if ii<6
        ylim([-.2 .6])
        plot([0 0],[-.2 .6],'--r');
        lgd = legend(name,'Location','northwest');
        lgd.FontSize = 8;
    else
        ylim([-.1 .2])
        plot([0 0],[-.1 .3],'--r');
    end
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


