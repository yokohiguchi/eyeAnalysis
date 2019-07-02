% Pupil analysis
% written by YH
% 6/11/2019

function D = eyeposition
clear; clear all;

% input
gType = input('Group type? (1:with training, 2:w/o training) :  ');
tType = input('Test type? (1:pre, 2:post, 3:re) :  ');
pl = input('Plot trial data? (y/n) :  ', 's');

if gType == 1
    ID = ['02';'03';'04';'05';'06';'08';'09';'11';'12';'14';'15';'16';'17';'21';'22';'24'];
    name = ['MH';'JY';'MN';'YM';'AU';'TM';'KK';'NW';'TW';'KT';'AS';'HK';'TI';'MK';'NK';'NF'];
    if tType == 1
        el_dir = 'data/eye_pre/';
        behav_dir = 'data/mat_pre/';
    elseif tType == 2
        el_dir = 'data/eye_post/';
        behav_dir = 'data/mat_post/';
    elseif tType == 3
        el_dir = 'data/eye_re/';
        behav_dir = 'data/mat_re/';
    end
    fig_dir = 'figure/eyeposi/training/';
else
    ID = ['01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16'];
    name = ['RC';'MT';'KO';'EI';'MY';'MM';'SM';'RU';'YU';'NW';'RF';'RK';'KM';'MM';'TM';'DS'];
    if tType == 1
        el_dir = 'data/eye_pre_ctrl/';
        behav_dir = 'data/mat_pre_ctrl/';
    elseif tType == 2
        el_dir = 'data/eye_post_ctrl/';
        behav_dir = 'data/mat_post_ctrl/';
    end
    fig_dir = 'figure/eyeposi/ctrl/';
end

nSubject = length(ID);

for sub = 1:nSubject
    
    fprintf('.');
    
    if pl=='y'
        pp = 1;
        h = figure(...
            'InvertHardcopy', 'off',...
            'Color', [1 1 1],...
            'Position', [0 0 1600 800]);
        hold all;
    end
    
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
                D.sub(sub).px(bb,nn)...
                    = mean(EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),6));
                D.sub(sub).py(bb,nn)...
                    = mean(EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),7));
                
                
                if pl=='y'
                    subplot(4,5,pp); pp=pp+1; hold all;
                    scatter(EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),6),...
                        EFIX(ismember(EFIX(:,1),find(cIdx==1 & eIdx==0)),7)); hold on;
                    scatter(mean(disk_x(cIdx==1 & eIdx==0)),mean(disk_y(cIdx==1 & eIdx==0)),300,'r');
                    scatter(Cfg.xCentre,Cfg.yCentre,100,'k+');
                    axis([0 Cfg.screensize_c 0 Cfg.screensize_r]);
                    set(gca,'YDir','reverse')
                    title(num2str(D.sub(sub).distance(bb,nn)))
                    drawnow;
                end
            else
                for train = [1 0]
                    nn=nn+1;
                    D.sub(sub).distance(bb,nn)...
                        = mean(EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),9));
                    D.sub(sub).px(bb,nn)...
                        = mean(EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),6));
                    D.sub(sub).py(bb,nn)...
                        = mean(EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),7));
                    
                    if pl=='y'
                        subplot(4,5,pp);  pp=pp+1; hold all;
                        scatter(EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),6),...
                            EFIX(ismember(EFIX(:,1),find(tIdx==train & cIdx==cond & eIdx==0)),7)); hold on;
                        scatter(mean(disk_x(tIdx==train & cIdx==cond & eIdx==0)),...
                            mean(disk_y(tIdx==train & cIdx==cond & eIdx==0)),300,'r');
                        scatter(Cfg.xCentre,Cfg.yCentre,100,'k+');
                        axis([0 Cfg.screensize_c 0 Cfg.screensize_r]);
                        set(gca,'YDir','reverse')
                        title(num2str(D.sub(sub).distance(bb,nn)))
                        drawnow;
                    end
                end
            end
        end
        
    end % for bb
    
    D.distance(sub,:) = mean(D.sub(sub).distance);
    D.px(sub,:) = mean(D.sub(sub).px);
    D.py(sub,:) = mean(D.sub(sub).py);
    
    if pl=='y'
        filename...
            = strcat(fig_dir,ID(sub,:),name(sub,:),'_',num2str(gType),'_',num2str(tType),'.fig');
        
        savefig(h,filename);
        close(h);
    end
    
end % for sub
px = D.px;
py = D.py;

px((mod(str2num(ID),4)+1==1),:)...
    =(Cfg.xCentre-px((mod(str2num(ID),4)+1==1),:))+Cfg.xCentre;
px((mod(str2num(ID),4)+1==3),:)...
    =(Cfg.xCentre-px((mod(str2num(ID),4)+1==3),:))+Cfg.xCentre;
py((mod(str2num(ID),4)+1==3),:)...
    =(Cfg.yCentre-py((mod(str2num(ID),4)+1==3),:))+Cfg.yCentre;
py((mod(str2num(ID),4)+1==4),:)...
    =(Cfg.yCentre-py((mod(str2num(ID),4)+1==4),:))+Cfg.yCentre;

diskLoc=[624,1296,624,1296;252,252,828,828];

figure(...
    'InvertHardcopy', 'off',...
    'Color', [1 1 1],...
    'Position', [0 0 1600 200]);
hold all;

for ii=1:5
    subplot(1,5,ii); hold all;
    scatter(px(:,ii),py(:,ii))
    scatter(mean(px(:,ii)),mean(py(:,ii)),'r','filled');
    if ii == 3 || ii == 5
        scatter(diskLoc(1,3),diskLoc(2,3),300,'r');
    else
        scatter(diskLoc(1,2),diskLoc(2,2),300,'r');
    end
    scatter(Cfg.xCentre,Cfg.yCentre,100,'k+');
    axis([0 Cfg.screensize_c 0 Cfg.screensize_r]);
    set(gca,'YDir','reverse')
    title(num2str(mean(D.distance(:,ii))))
end
return