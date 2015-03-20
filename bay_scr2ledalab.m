function out = bay_scr2ledalab(whattodo,varargin)

% whattodo: 
%   'run' - runs Ledalab analysis for chosen sessions & subjects CAVE:
%           define Ledalab settings before
%   'load_results' - loads results data
%   'plot' - produces .png of results data (phasic, tonic, driver, impulse 
%            response funct., condition information)
%   'get_timecourse' - produces binning matrix per subject 3 dimensional
%                      (defined timepoints (target_timebin), trials, 
%                      amplitude) -> amplitudes cluster matrix per bin & trial

addpath('/home/grahl/Scripts/globalfunctions/trunk');
addpath(genpath('/home/grahl/Scripts/Ledalab/trunk'));
%data path
path         =  '/common/prd03/projects/bay/behavior/BayBP/SCR/toLedalab/';
out = [];
%% read and make Ledalab analysis
if strcmp(whattodo,'run');
 
    sessions = varargin{1};
    subjects = varargin{2};
    
    for session   = sessions;
        for subject  = subjects
            fprintf('Running: Session: %03d; Subject: %03d...\n',session,subject);
            try
                fname        = sprintf('%stoLedalab_BayBP_S%d_%02d.mat',path,session,subject);
                Ledalab({fname}, 'open', 'mat','analyze','CDA', 'optimize',6, 'overview',  1, 'export_era', [0 18 0 1], 'export_scrlist', [0 1], 'export_eta', 1);
            end
        end
    end

%% load results data
elseif strcmp(whattodo,'load_results');
    %
    session      = varargin{1};
    subject      = varargin{2};
    fname        = sprintf('%stoLedalab_BayBP_S%d_%02d_results.mat',path,session,subject);
    if exist(fname) ~= 0
        out          = load(fname);
    else
        fprintf('File %s doesn''t exist.\n',fname);
        return
    end

%% plot SCR analysis (Ledalab)
elseif strcmp(whattodo,'plot');
    close all
    X = 40;                    % x-axis paper size
    Y = 21;                    % y-axis paper size
    xMargin = 2;               % left/right margins from page borders
    yMargin = 1;               % bottom/top margins from page borders
    xSize = X - 2*xMargin;     % figure size on paper (width & height)
    ySize = Y - 2*yMargin;     % figure size on paper (width & height)
    F = figure;
    set(F, 'PaperUnits','centimeters')
    set(F, 'PaperSize',[X Y])
    set(F, 'PaperPosition',[xMargin yMargin xSize ySize])

    subplot('position',[0.025 0.15 0.95 0.8])
    session      = varargin{1};
    subject      = varargin{2};
    d            = bay_scr2ledalab('load_results',session,subject);
    fname        = sprintf('%stoLedalab_BayBP_S%d_%02d_plot.png',path,session,subject);
    fprintf('Running: Session: %03d; Subject: %03d...\n',session,subject);
    if ~isempty(d)
        %plot the tonic and phasic drive
        plot(d.data.time.data, d.data.conductance.data,'k.-');%the raw data
        hold on
        plot(d.data.time.data, d.analysis.phasicData,'r','linewidth',2);
        plot(d.data.time.data, d.analysis.tonicData ,'m','linewidth',2);
        plot(d.data.time.data(1:length( d.analysis.kernel)),d.analysis.kernel);
        m = d.analysis.driver;
        %plot(data.time,zscore(analysis.driver)*std(analysis.phasicData)+m,'b','linewidth',1);
        h = area(d.data.time.data,zscore(d.analysis.driver)*std(d.analysis.phasicData)+m,0,'facecolor','b','linewidth',1,'edgecolor','b');
        alpha(0.3)
        box off
        %
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        axis tight
        title(sprintf('Subject: %d, Phase: %d, CompoundError: %0.3g,Negativity:%0.3g',subject,session,d.analysis.error.negativity,d.analysis.error.compound))
        
        %% plot a line at stimulus onset.
        tface = 4;
        colors = [0 0 0; 1 0 0; .5 .5 .5; 0 0 1];
        for ne = 1:length(d.data.event);
            %find the color that corresponds to the current contiion
            cond = d.data.event(ne).nid;
            if isnan(cond) == 0
                color = colors( cond,:);
            else
                color = [ 0 0 0];
                cond  = 10;
            end
            %plot stimulus onsets
            patchline(repmat(d.data.event(ne).time - d.data.time.timeoff,1,2),ylim,'edgecolor',color,'linewidth',2,'edgealpha',0.25);
        end
        xlabel('time (seconds)')
        
        %% plot the kernel as an inset
        hold on
        axes('position',[0.75 0.75 0.2 0.2])
        
        %plot the bateman (tau 1 & 2)
        plot(d.data.time.data(1:length( d.analysis.kernel)),d.analysis.kernel,'k','linewidth',2);
        title(sprintf('Bateman(%0.3g,%0.3g)',d.analysis.tau(1),d.analysis.tau(2)))
        xlabel('time (s)')
        box off
        
        %% save the stuff
        SaveFigure(fname) % adapt saving format in SaveFigure (jpg, pdf, png,...)
    end

%% build 3D matrix -> SCR amplitude, time bins, trials (here: 22 heat onsets)
elseif strcmp(whattodo,'get_timecourse')
    session        = varargin{1};
    target_timebin = varargin{3};
    counter        = 0;
    Mbinned        = [];
    
    for subject    = varargin{2};
        fprintf('Processing subject #%03d...\n',subject);
        counter = counter + 1;
        %extract the data triggered with heat onset
        out     = bay_scr2ledalab('load_results',session,subject);
        if ~isempty(out)
            M_trial       = out.analysis.split_driver.y; % deconvolved SCR
            c       = out.analysis.split_driver.c; % trigger information
            M_trial       = M_trial(:,ismember(c,[2 4]));
            %create a binning matrix
            [BM,centers]=BinningMatrix(size(M_trial,1),target_timebin);
            Mbinned(:,:,counter)       = BM'*M_trial;
        end
    end
    fname = sprintf('%stoLedalab_BayBP_S%d_%02d_matrix.mat',path,session,subject);
    out = Mbinned;
    save(fname, 'out');

%% plot all subjects heat matrices within one figure (per session)
elseif strcmp(whattodo,'plot_matrix') % insert session, subjects, time-bin value, time-window in sec
    session         = varargin{1};
    subject         = varargin{2};
    target_timebin  = varargin{3};
    time_window_sec = varargin{4}; % extracted time-window in sec. after trigger onset
    plotSizeVer     = 4; % subplot size vertical
    plotSizeHor     = 5; % subplot size horizontal
    
    close all
    X = 50;                  % x-axis paper size
    Y = 25;                    % y-axis paper size
    xMargin = .1;               % left/right margins from page borders
    yMargin = .2;               % bottom/top margins from page borders
    xSize = X - 2*xMargin;     % figure size on paper (width & height)
    ySize = Y - 2*yMargin;     % figure size on paper (width & height)
    F = figure;
        
    set(F, 'PaperUnits','centimeters');
    set(F, 'PaperSize',[X Y]);
    set(F, 'PaperPosition',[xMargin yMargin xSize ySize]);
    
    subplot(plotSizeVer,plotSizeHor,1);
    set(gca,'xtick',[],'ytick',[]);
    axis off;
    % adapt information text
        text(-.2,1.05,sprintf('SCR amplitudes - %02d sec. heat onset\n1 time-bin: %03d ms\ntrials: 1-11 no TENS, 12-22 TENS\nheat plateau after approx. 1 sec.',time_window_sec,(target_timebin/time_window_sec)*10), ...
        'Units', 'normalized', ... 
        'VerticalAlignment', 'top', ... 
        'HorizontalAlignment', 'left', ...
        'Color', 'black',...
        'FontWeight','bold',...
        'FontSize',14);
    
    s = 1;
    for subject = varargin{2};
        s = s+1;
        load(sprintf('%stoLedalab_BayBP_S%d_%02d_matrix.mat',path,session,subject));
        % bay_scr2ledalab('get_timecourse',session,subject,target_timebin);
        fprintf('Processing subject #%03d...\n',subject);
        subplot(plotSizeVer,plotSizeHor,s);
        imagesc(out);
%         imagesc(out(1:8,:)); % select time-bins you want to dislay in matrices
        colorbar;
        xlabel('trials','FontWeight','bold');
        ylabel('time','FontWeight','bold');
        title(sprintf('BayBP_S%d_%02d',session,subject),'interpreter','none','FontWeight','bold','FontSize', 16);
        
    end
    fname = sprintf('%sBayBP_S%d_heatMatrices_%03dbins_heatOn%02dsec.pdf',path,session,target_timebin,time_window_sec);
    SaveFigure(fname); % adapt saving format in SaveFigure (jpg, pdf, png,...)
        
    
    
%% compute mean and variance  
elseif strcmp(whattodo,'calc_stats')
    session         = varargin{1};
    subject         = varargin{2};
    for session = varargin{1}
        for subject = varargin{2};
            load(sprintf('%stoLedalab_BayBP_S%d_%02d_matrix.mat',path,session,subject));
            fprintf('Processing subject #%03d...\n',subject);
            % calculate
            M_trial = mean(out); % per trial across bins (1 value per column)
            Variance_trial = var(out);
            SD_trial = std(out);
            M_time = mean(out,2); % per bin across trials (1 value per row)
            Variance_time = var(out,0,2);
            SD_time = std(out,0,2);
            % save
            s.mean_trial = M_trial;
            s.variance_trial = Variance_trial;
            s.standDev_trial = SD_trial;
            s.mean_time = M_time;
            s.variance_time = Variance_time;
            s.standDev_time = SD_time;
        end
        if session == 1 S1 = [s]; end;
        if session == 2 S2 = [s]; end;
        
    end
    
    fname = sprintf('%stoLedalab_BayBP_%02d_descrStats.mat',path,subject);
    descr_stats.S1 = S1;
    descr_stats.S2 = S2;
    save(fname, 'descr_stats');
    
end

% % bay_scr2ledalab('run',1:2,[1:3 5:19]);
% % clear all; close all;
% % bay_scr2ledalab('load_results',1:2,[1:3 5:19]);
% % clear all; close all;
% % bay_scr2ledalab('plot',1:2,[1:3 5:19]);
% % clear all; close all;
% % bay_scr2ledalab('get_timecourse',1:2,[1:3 5:19]);
% % clear all; close all;


% % %%
% % clear all;
% % path         =  '/Users/onat/Desktop/toLedalab/corrected/';
% % session = 1;
% % for subject = 2:19;
% %     try
% %         fname   = sprintf('%stoLedalab_BayBP_S%d_%02d_results.mat',path,session,subject);
% %         load(fname);
% %         d = analysis.split_driver
% %         figure;
% %         plot(d.x(:,1:4),d.mean)
% %     end
% % end
% % %%
% %
% % %%
% % %color code for tens vs. no-tens
% % co{[43]} = 'r';
% % co{[44]} = 'b';
% % co{[23]} = 'r';
% % co{[24]} = 'b';
% % %
% % t{[1]}   = { 'cue'};
% % t{[21]}  = { 'heat'};
% % c        = 0;
% % %%
% % figure(100);
% % c    = 0;
% % for cue = [ 0 20];
% %     c = c + 1;
% %     subplot(1,2,c);
% %     for block = [23 24]+cue
% %         i = analysis.split_phasicData.c == block ;
% %         plot(mean(analysis.split_phasicData.y(:,i),2),'color',co{block})
% %         hold on;
% %     end
% %     [cue block]
% %     title(t{cue+1});
% %     legend('tens-','tens+');
% %     hold off
% % end
