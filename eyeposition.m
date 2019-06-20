% Pupil analysis
% written by YH
% 6/11/2019

function D = eyeposition
clear; clear all;

% Load data
load('02MH13yh_eye.mat')
load('02_MH_1_3.mat')

nn = 0;
for block = 1:length([Data.condition])
    for trl = 1:length(Data(block).TR)
        nn = nn + 1;
        tIdx(nn) = Data(block).TR(trl).train;
        cIdx(nn) = Data(block).condition;
        eIdx(nn) = Data(block).TR(trl).eyeError;
        disk_x(nn) = Data(block).TR(trl).dispLocPeriph(1);
        disk_y(nn) = Data(block).TR(trl).dispLocPeriph(2);
    end
end

% EFIX(trl,0,stime,etime,dur,x,y,pupil,dist)
for ii=1:length(EFIX)
    EFIX(ii,9)...
        = (sqrt((EFIX(ii,6) - disk_x(EFIX(ii,1)))^2 ...
        + (EFIX(ii,7) - disk_y(EFIX(ii,1)))^2))/Cfg.pixelsPerDegree;
end

figure(...
    'InvertHardcopy', 'off',...
    'Color', [1 1 1],...
    'Position', [0 0 1600 200]);
nn=0;
for cond = 1:3
    if cond == 1
        nn=nn+1;
        D.distance(nn)...
            = mean(EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),9));
        subplot(1,5,nn);
        scatter(EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),6),...
            EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),7)); hold on;
        scatter(mean(disk_x(cIdx==1 & eIdx==0)),mean(disk_y(cIdx==1 & eIdx==0)),300,'r');
        axis([0 Cfg.screensize_c 0 Cfg.screensize_r]);
        set(gca,'YDir','reverse')
        title(num2str(D.distance(nn)))
    else
        for train = [1 0]
            nn=nn+1;
            D.distance(nn)...
            = mean(EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),9));
        subplot(1,5,nn);
        scatter(EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),6),...
            EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),7)); hold on;
        scatter(mean(disk_x(tIdx==train & cIdx==cond & eIdx==0)),...
            mean(disk_y(tIdx==train & cIdx==cond & eIdx==0)),300,'r');
        axis([0 Cfg.screensize_c 0 Cfg.screensize_r]);
        set(gca,'YDir','reverse')
        title(num2str(D.distance(nn)))
        end
    end
end


return