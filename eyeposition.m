% Pupil analysis
% written by YH
% 6/11/2019

function D = eyeposition_pre
clear; clear all;


ID = ['02';'03';'04';'05';'06';'08';'09';'11';'12';'14';'15';'16';'17';'21';'22';'24'];
name = ['MH';'JY';'MN';'YM';'AU';'TM';'KK';'NW';'TW';'KT';'AS';'HK';'TI';'MK';'NK';'NF'];

nSubject = length(ID);

for sub = 1:nSubject
    
    pp = 1;
    figure(...
        'InvertHardcopy', 'off',...
        'Color', [1 1 1],...
        'Position', [0 0 1600 800]);
    hold all;
    
    for bb = 1:4
        if bb == 1
            sess = 1; r = 3;
        else
            sess = 2; r = bb-1;
        end
        
        % Load data
        el...
            = strcat('data/eye_pre/',ID(sub,:),name(sub,:),...
            num2str(sess),num2str(r),'yh_eye.mat');
        behav...
            = strcat('data/mat_pre/',ID(sub,:),'_',name(sub,:),...
            '_',num2str(sess),'_',num2str(r),'.mat');
        
        load(el);
        load(behav);
        
        % Index
        ii=0;
        for block = 1:length([Data.condition])
            for trl = 1:length(Data(block).TR)
                ii = ii + 1;
                tIdx(ii) = Data(block).TR(trl).train;
                cIdx(ii) = Data(block).condition;
                eIdx(ii) = Data(block).TR(trl).eyeError;
                disk_x(ii) = round(Data(block).TR(trl).dispLocPeriph(1));
                disk_y(ii) = round(Data(block).TR(trl).dispLocPeriph(2));
            end
        end
        clear ii
        
        % EFIX(trl,0,stime,etime,dur,x,y,pupil,dist)
        % *Distance from center to disk: 9.22
        for ii=1:length(EFIX)
            EFIX(ii,9)...
                = (sqrt((EFIX(ii,6) - disk_x(EFIX(ii,1)))^2 ...
                + (EFIX(ii,7) - disk_y(EFIX(ii,1)))^2))/Cfg.pixelsPerDegree;
            
            %             % Find fixations after disk mask
            %             if EFIX(ii,3)<= MSGtime((MSGtime(:,1)==EFIX(ii,1) & MSGtime(:,3)==20),2)
            %                 EFIX(ii,10)=1;
            %             else
            %                 EFIX(ii,10)=0;
            %             end
            
%             % Find fixations 3 degrees away from the center
%             EFIX(ii,10)...
%                 = (sqrt((EFIX(ii,6) - Cfg.xCentre)^2 ...
%                 + (EFIX(ii,7) - Cfg.yCentre)^2))/Cfg.pixelsPerDegree;
            
        end
        clear ii
        
%         % Remove fixations 3 degrees away from the center
%         EFIX((EFIX(:,10)>3),:)=[];
        
        nn=0;
        for cond = 1:3
            if cond == 1
                nn=nn+1;
                D.sub(sub).distance(bb,nn)...
                    = mean(EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),9));
                subplot(4,5,pp); pp=pp+1; hold all;
                scatter(EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),6),...
                    EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),7)); hold on;
                scatter(mean(disk_x(cIdx==1 & eIdx==0)),mean(disk_y(cIdx==1 & eIdx==0)),300,'r');
                scatter(Cfg.xCentre,Cfg.yCentre,100,'k+');
                axis([0 Cfg.screensize_c 0 Cfg.screensize_r]);
                set(gca,'YDir','reverse')
                title(num2str(D.sub(sub).distance(bb,nn)))
            else
                for train = [1 0]
                    nn=nn+1;
                    D.sub(sub).distance(bb,nn)...
                        = mean(EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),9));
                    subplot(4,5,pp);  pp=pp+1; hold all;
                    scatter(EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),6),...
                        EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),7)); hold on;
                    scatter(mean(disk_x(tIdx==train & cIdx==cond & eIdx==0)),...
                        mean(disk_y(tIdx==train & cIdx==cond & eIdx==0)),300,'r');
                    scatter(Cfg.xCentre,Cfg.yCentre,100,'k+');
                    axis([0 Cfg.screensize_c 0 Cfg.screensize_r]);
                    set(gca,'YDir','reverse')
                    title(num2str(D.sub(sub).distance(bb,nn)))
                end
            end
        end
        
    end % for bb
    
    D.distance(sub,:) = mean(D.sub(sub).distance);
    
end % for sub
return