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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ledalab functions:
% open: biotrace, biopac, biopacmat, cassylab, leda (default), varioport, visionanalyzer, vitaport, portilab, psychlab, mat, text (Type 1), text2 (Type 2)
% filter: [filter_order, lower-cutoff frequency]; default = no filter applied
% downsample: 0 (default) , 2+
% smooth: {type, width} -> type: mean (moving average), hann (hanning window), gauss (gauss window), adapt (adaptive smoothing using gauss); 
%                       -> width: width of smoothing window in samples (does not apply for adapt)
% analyze: 'CDA' (Discrete Decomposition Analysis), 'DDA' (Continuous Decomposition Analysis)
% optimize: 2 (default), 1 - 6
% export_era: [respwin_start respwin_end amp_threshold (filetype)]; filetype: 1 (Matlab file = default), 2 (Text file), 3 (Excel file)
% export_scrlist: [amp_threshold (filetype)]; filetype: 1 (Matlab file = default), 2 (Text file), 3 (Excel file)
% overview: 0 (default), 1 -> Save decomposition result to jpg for easy control (boolean option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; close all;
% bay_scr2ledalab('run_one',1:2,[1:3 5:19]);
% clear all; close all;
% bay_scr2ledalab('run_two',20:35);
% clear all; close all;
% bay_scr2ledalab('run_one',1:2,4);
% clear all; close all;


addpath('/home/grahl/Scripts/globalfunctions/trunk');
addpath(genpath('/home/grahl/Scripts/Ledalab/trunk'));
%data path
path         =  '/common/prd03/projects/bay/behavior/BayBP/SCR/toLedalab/';
pathHome     =  '/home/grahl/Bay/Exp_bay/BehaviorP/Pics/';
out = [];
%% read and make Ledalab analysis with one session per file
if strcmp(whattodo,'run_one');
 
    sessions = varargin{1};
    subjects = varargin{2};
    
        for session  = sessions;
        for subject  = subjects
            fprintf('Running: Subject: %03d... Session: %d...\n',subject,session);
            try
                fname        = sprintf('%stoLedalab_BayBP_S%d_%02d.mat',path,session,subject);
%                 fname        = sprintf('%stoLedalab_BayBP_%02d.mat',path,subject);
                Ledalab({fname}, 'open', 'mat','analyze','CDA', 'optimize',6, 'overview',  1, 'export_era', [0 12 0 1], 'export_scrlist', [0 1], 'export_eta', 1);
            end
        end
        end
        
%% read and make Ledalab analysis with two sessions per file
elseif strcmp(whattodo,'run_two');
    path     =  '/common/prd03/projects/bay/behavior/BayBP/SCR/toLedalab/BayBP_20to35mat/';
    subjects = varargin{1};
    
        for subject  = subjects
            fprintf('Running: Subject: %03d...\n',subject);
            try
                fname        = sprintf('%stoLedalab_BayBP_%02d.mat',path,subject);
                Ledalab({fname}, 'open', 'mat','analyze','CDA', 'optimize',6, 'overview',  1, 'export_era', [0 12 0 1], 'export_scrlist', [0 1], 'export_eta', 1);
            end
        end

%% load results data
elseif strcmp(whattodo,'load_results');
    %
    session      = varargin{1};
    subject      = varargin{2};
    if session   == 0
        fname    = sprintf('%stoLedalab_BayBP_%02d_results.mat',path,subject);
    else
        fname    = sprintf('%stoLedalab_BayBP_S%d_%02d_results.mat',path,session,subject);
    end
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
    if session   == 0
        d        = bay_scr2ledalab('load_results',session,subject);
    else
        d        = bay_scr2ledalab('load_results',session,subject);
    end
    if session   == 0
        fname    = sprintf('%stoLedalab_BayBP_%02d_plot',path,subject); 
    else
        fname    = sprintf('%stoLedalab_BayBP_S%d_%02d_plot',path,session,subject);
    end
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
        colors = [0 0 0; 1 0 0; .5 .5 .5; 0 0 1; 0 1 0];
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
%     1 - cue without TENS
%     2 - heat without TENS
%     3 - cue with TENS
%     4 - heat with TENS
    session        = varargin{1};
    target_timebin = varargin{3};
    counter        = 0;
    Mbinned        = [];
    
        for subject    = varargin{2};
            fprintf('Processing session #%1d subject #%03d...\n',session,subject);
            counter = 1; % saves one file per subject
            %counter + 1; % builds N dimensional matrix (N = subject size)
            clear out;
            %extract the data triggered with heat onset
            if subject < 20
                out = bay_scr2ledalab('load_results',session,subject);
            elseif subject > 19
                out = bay_scr2ledalab('load_results',0,subject);
            end
            if ~isempty(out)
                if subject < 20
                    M_trial = out.analysis.split_driver.y; % deconvolved SCR
                    c       = out.analysis.split_driver.c; % trigger information
                    M_trial = M_trial(:,ismember(c,[2 4]));
                elseif subject > 19
                    if session == 1
                        M_trial = out.analysis.split_driver.y(:,[3:13 25:35 47:57 69:79]); % deconvolved SCR columns of session 1
                        c       = out.analysis.split_driver.c(:,[3:13 25:35 47:57 69:79]); % trigger information columns of session 1
                        M_trial = M_trial(:,ismember(c,[2 4]));
                    elseif session == 2
                        M_trial = out.analysis.split_driver.y(:,[14:24 36:46 58:68 80:90]); % deconvolved SCR columns of session 2
                        c       = out.analysis.split_driver.c(:,[14:24 36:46 58:68 80:90]); % trigger information columns of session 2
                        M_trial = M_trial(:,ismember(c,[2 4]));
                    end
                end
                %create a binning matrix
                [BM,centers] = BinningMatrix(size(M_trial,1),target_timebin);
                
                if subject == 4 && session == 1
                    Mbinned(:,2:21,counter) = BM'*M_trial;
                    Mbinned(:,1,counter) = NaN;
                    Mbinned(:,22,counter) = NaN;
                else
                    Mbinned(:,:,counter) = BM'*M_trial;
                end
            end
            fname = sprintf('%stoLedalab_BayBP_S%d_%02d_matrix.mat',path,session,subject);
            out = Mbinned;
            save(fname, 'out');
        end

%% plot all subjects heat matrices within one figure (per session)
elseif strcmp(whattodo,'plot_matrix') % insert session, subjects, time-bin value, time-window in sec
    session         = varargin{1};
    subject         = varargin{2};
    target_timebin  = varargin{3};
    time_window_sec = varargin{4}; % extracted time-window in sec. after trigger onset
    plotSizeVer     = 7; % subplot size vertical
    plotSizeHor     = 5; % subplot size horizontal
    
    close all
    X = 70;                  % x-axis paper size
    Y = 40;                    % y-axis paper size
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
    fname = sprintf('%sBayBP_S%d_heatMatrices_%03dbins_heatOn%02dsec_%02dsubs.pdf',path,session,target_timebin,time_window_sec,s-1);
    SaveFigure(fname); % adapt saving format in SaveFigure (jpg, pdf, png,...)
        
    
    
%% compute mean and variance  
elseif strcmp(whattodo,'calc_stats')
    session         = varargin{1};
    subject         = varargin{2};
    target_timebin  = varargin{3};
    time_window_sec = varargin{4};
    for subject = varargin{2}
        for session = varargin{1};
            load(sprintf('%stoLedalab_BayBP_S%d_%02d_matrix.mat',path,session,subject));
            fprintf('Processing subject #%03d... session #%d...\n',subject,session);
            % calculate descriptive statistics
            % per trial across bins/time
            M_trial = nanmean(out); % 1 value per column/trial
            Variance_trial = nanvar(out);
            SD_trial = nanstd(out);
            % per bin across 11 trials per condition
            % first column = trigger 2 (heat no TENS), second column = trigger 4 (heat TENS)
            M_time(:,1) = nanmean(out(:,1:11),2); % 1 value per row per condition -> 2 columns with n bin rows
            M_time(:,2) = nanmean(out(:,12:22),2);
            Variance_time(:,1) = nanvar(out(:,1:11),0,2); % variance across trials
            Variance_time(:,2) = nanvar(out(:,12:22),0,2);
            SD_time(:,1) = nanstd(out(:,1:11),0,2); % SD across trials
            SD_time(:,2) = nanstd(out(:,12:22),0,2);
            % statistics_trial: Test the null hypothesis that the variances per trial across 
            % n time-bins are  equal comparing heatNoTens-trial with corresponding heatTens- 
            % trial, using the Brown-Forsythe test (uses median instead of mean!)
            for k = 1:11
                [p_trial(k,1),stats_trial(k,1)] = vartestn(out(:,[k k+11]),'TestType','BrownForsythe','Display','off');
%                 [p_trial(k,1),stats_trial(k,1)] = vartest(out(:,[k k+11]));
                statsF_trial(k,1) = stats_trial(k).fstat(:);
                [~,p2_trial(k,1:2),CI(k,1:2),statsT1_trial(k,1)] = ttest2(out(:,k),out(:,k+11));
                statsT2_trial(k,1) = statsT1_trial(k).tstat;
%                 [~,p2_trial(k,1:2),CI(k,1:2),statsT1_trial(k,1)]=ttest(out(:,k),out(:,k+11));
            end
            % statistics_time: Test the null hypothesis that the variances per time-bin across 
            % n trials are  equal comparing heatNoTens-time-bin with corresponding heatTens- 
            % time-bin, using the Brown-Forsythe test (uses median instead of mean!)
            out2 = out';
            for k = 1:target_timebin
                clear temp_time;
                temp_time = [out2(1:11,k) out2(12:22,k)];
                [p_time(k,1),stats_time(k,1)] = vartestn(temp_time,'TestType','BrownForsythe','Display','off');
%                 [p_time(k,1),stats_time(k,1)] = vartest(temp_time);
                statsF_time(k,1) = stats_time(k).fstat(:);
                
            end
            % save
            s.mean_trial = M_trial;
            s.variance_trial = Variance_trial;
            s.standDev_trial = SD_trial;
            s.mean_time = M_time;
            s.variance_time = Variance_time;
            s.standDev_time = SD_time;
            s.ftestPval_trial = p_trial;
            s.ftestF_trial = stats_trial;
            s.ftestFval_trial = statsF_trial;
            s.ftestPval_time = p_time;
            s.ftestF_time = stats_time;
            s.ftestFval_time = statsF_time;
            
            s.ttestPval_trial = p2_trial;
            s.ttestT_trial = statsT1_trial;
            s.ttestTval_trial = statsT2_trial;
            s.ttestPval_time = p2_time;
            s.ttestT_time = statsT_time;
            s.ttestTval_time = statsT_time;
            
            if session == 1  descr_stats.S1 = [s]; end;
            if session == 2  descr_stats.S2 = [s]; end;
            
        end
        fname = sprintf('%stoLedalab_BayBP_%02d_descrStats.mat',path,subject);
        save(fname, 'descr_stats');
        
    end
    
    
    
    
    
%% compute mean and variance across subjects & plot figure
elseif strcmp(whattodo,'calc_allSubsStats')
    target_timebin = varargin{2};
    % session        = varargin{1};
    counter        = 0;
    statsAll       = [];
    for subject    = varargin{1};
        fprintf('Processing subject #%03d...\n',subject);
        counter = counter + 1;
        %extract the data triggered with heat onset
        load(sprintf('%stoLedalab_BayBP_%02d_descrStats.mat',path,subject));
        %create a binning matrix
        statsAll.S1.mean_trial(:,:,counter) = descr_stats.S1.mean_trial;
        statsAll.S1.mean_time(:,:,counter)  = descr_stats.S1.mean_time';
        statsAll.S2.mean_trial(:,:,counter) = descr_stats.S2.mean_trial;
        statsAll.S2.mean_time(:,:,counter)  = descr_stats.S2.mean_time';
        statsAll.S1.SD_trial(:,:,counter) = descr_stats.S1.standDev_trial;
        statsAll.S1.SD_time(:,:,counter)  = descr_stats.S1.standDev_time';
        statsAll.S2.SD_trial(:,:,counter) = descr_stats.S2.standDev_trial;
        statsAll.S2.SD_time(:,:,counter)  = descr_stats.S2.standDev_time';
        statsAll.S1.variance_trial(:,:,counter) = descr_stats.S1.variance_trial;
        statsAll.S1.variance_time(:,:,counter)  = descr_stats.S1.variance_time';
        statsAll.S2.variance_trial(:,:,counter) = descr_stats.S2.variance_trial;
        statsAll.S2.variance_time(:,:,counter)  = descr_stats.S2.variance_time';
        statsAll.S1.ftestFval_trial(:,:,counter)= descr_stats.S1.ftestFval_trial;
        statsAll.S1.ftestFval_time(:,:,counter) = descr_stats.S1.ftestFval_time;
        statsAll.S2.ftestFval_trial(:,:,counter)= descr_stats.S2.ftestFval_trial;
        statsAll.S2.ftestFval_time(:,:,counter) = descr_stats.S2.ftestFval_time;
        
    end

    fname = sprintf('%sBayBP_descrStatsAll.mat',path);
    save(fname, 'statsAll');
    
    % plot all
    % predefine settings
    subs                = 19;
    subs_con = [3 5 7 12 14 16 19]; % [1 3 5 7 10 12 14 16 17 19];
    subs_var = [2 4 6 8 9 11 13 15 18]; % [2 4 6 8 9 11 13 15 18];
    x_range_trial       = [0 12];
    x_range_time        = [-10 190];
    y_range_trial       = [0 1050];
    y_range_time        = [0 1200];
    xLabel_plots_trial  = 'trials';
    xLabel_plots_time   = 'time (in 100 ms)';
    yLabel_plots        = 'SCR mean amplitude in \muS (+/- STD)';
    legend_label        = {'no TENS', 'TENS'};
    linePlot            = 2; % linewidth of lines within plots
    colConN             = [1 0 0]; % red
    colConT             = [1 0.6 0.6]; % light red
    colVarN             = [0 0 1]; % blue
    colVarT             = [0.6 0.6 1]; % light blue
    X = 40;                  % A4 paper size
    Y = 21;                    % A4 paper size
    xMargin = 1;               % left/right margins from page borders
    yMargin = 1;               % bottom/top margins from page borders
    xSize = X - 2*xMargin;     % figure size on paper (width & height)
    ySize = Y - 2*yMargin;     % figure size on paper (width & height)
    
    % each subplot = 2 figures
    % define axis properties
    width_fig1 = 0.3;
    height_fig = 0.35;
    width_fig2 = 0.10;
    X_left1 = 0.05;
    X_left2 = 0.38;
    X_right1 = 0.54;
    X_right2 = 0.87;
    Y_top = 0.58;
    Y_bottom = 0.06;
    Line_with = 2;
    y = 0.99; % textbox coordinates
    x = 0.35;
    
    % calculate mean of all F-test F-values comparing heat-NoTens vs. heat-Tens 
    % -> per trial:
    F_S1_trial_con = nanmean(statsAll.S1.ftestFval_trial(:,:,subs_con),3);
    F_S1_trial_var = nanmean(statsAll.S1.ftestFval_trial(:,:,subs_var),3);
    F_S2_trial_con = nanmean(statsAll.S2.ftestFval_trial(:,:,subs_con),3);
    F_S2_trial_var = nanmean(statsAll.S2.ftestFval_trial(:,:,subs_var),3);
     % -> per time-bin:
    F_S1_time_con = nanmean(statsAll.S1.ftestFval_time(:,:,subs_con),3);
    F_S1_time_var = nanmean(statsAll.S1.ftestFval_time(:,:,subs_var),3);
    F_S2_time_con = nanmean(statsAll.S2.ftestFval_time(:,:,subs_con),3);
    F_S2_time_var = nanmean(statsAll.S2.ftestFval_time(:,:,subs_var),3);
    
    % plot f-stats of Brown-Forsythe (2x2 subplots)
    f0 = figure;
    set(f0, 'PaperUnits','centimeters');
    set(f0, 'PaperSize',[X Y]);
    set(f0, 'PaperPosition',[xMargin yMargin xSize ySize]);
    set(f0, 'PaperOrientation','Portrait');
    legend_f0 = {sprintf('constant (N = %.0f)',length(subs_con)),sprintf('variable (N = %.0f)',length(subs_var))};
    subplot(2,2,1);plot(F_S1_trial_con,'r');hold on; plot(F_S1_trial_var,'b');
    legend(legend_f0); legend boxoff;
    xlabel(xLabel_plots_trial,'FontWeight','bold','FontSize',12);
    ylabel('mean F-value(1,358)','FontWeight','bold','FontSize',12);
    title('S1, Brown-Forsythe (homoscedasticity): Tens vs. noTens','FontWeight','bold','FontSize',14);
    xlim(x_range_trial); ylim([0 250]);hold off;
    subplot(2,2,2);plot(F_S2_trial_con,'r');hold on; plot(F_S2_trial_var,'b');hold off;
    legend(legend_f0); legend boxoff;
    xlabel(xLabel_plots_trial,'FontWeight','bold','FontSize',12);
    ylabel('mean F-value(1,358)','FontWeight','bold','FontSize',12);
    title('S2, Brown-Forsythe (homoscedasticity): Tens vs. noTens','FontWeight','bold','FontSize',14);
    xlim(x_range_trial); ylim([0 250]);hold off;
    subplot(2,2,3);plot(F_S1_time_con,'r');hold on; plot(F_S1_time_var,'b');
    legend(legend_f0); legend boxoff;
    xlabel(xLabel_plots_time,'FontWeight','bold','FontSize',12);
    ylabel('mean F-value(1,20)','FontWeight','bold','FontSize',12);
    title('S1, Brown-Forsythe (homoscedasticity): Tens vs. noTens','FontWeight','bold','FontSize',14);
    xlim(x_range_time); ylim([0 6]);hold off;
    subplot(2,2,4);plot(F_S2_time_con,'r');hold on; plot(F_S2_time_var,'b');hold off;
    legend(legend_f0); legend boxoff;
    xlabel(xLabel_plots_time,'FontWeight','bold','FontSize',12);
    ylabel('mean F-value(1,20)','FontWeight','bold','FontSize',12);
    title('S2, Brown-Forsythe (homoscedasticity): Tens vs. noTens','FontWeight','bold','FontSize',14);
    xlim(x_range_time); ylim([0 6]);hold off;
    
    
    
    f1 = figure;
    set(f1, 'PaperUnits','centimeters');
    set(f1, 'PaperSize',[X Y]);
    set(f1, 'PaperPosition',[xMargin yMargin xSize ySize]);
    set(f1, 'PaperOrientation','Portrait');
    
    % subplot 1 - constant conditioning
    % [startX startX width height]
    a1 = axes('position', [X_left1 Y_top width_fig1 height_fig]);
    a2 = axes('position', [X_left2 Y_top width_fig2 height_fig]);
    a3 = axes('position', [X_right1 Y_top width_fig1 height_fig]);
    a4 = axes('position', [X_right2 Y_top width_fig2 height_fig]);
    % subplot 3 - variable conditioning
    a5 = axes('position', [X_left1 Y_bottom width_fig1 height_fig]);
    a6 = axes('position', [X_left2 Y_bottom width_fig2 height_fig]);
    % subplot 4 - variable test phase
    a7 = axes('position', [X_right1 Y_bottom width_fig1 height_fig]);
    a8 = axes('position', [X_right2 Y_bottom width_fig2 height_fig]);
    
    % calculate mean -> 1 value per variable
    M_S1_trial_conN = nanmean(nanmean(statsAll.S1.mean_trial(:,1:11,subs_con)));
    M_S1_trial_conT = nanmean(nanmean(statsAll.S1.mean_trial(:,12:22,subs_con)));
    M_S2_trial_conN = nanmean(nanmean(statsAll.S2.mean_trial(:,1:11,subs_con)));
    M_S2_trial_conT = nanmean(nanmean(statsAll.S2.mean_trial(:,12:22,subs_con)));
    M_S1_trial_varN = nanmean(nanmean(statsAll.S1.mean_trial(:,1:11,subs_var)));
    M_S1_trial_varT = nanmean(nanmean(statsAll.S1.mean_trial(:,12:22,subs_var)));
    M_S2_trial_varN = nanmean(nanmean(statsAll.S2.mean_trial(:,1:11,subs_var)));
    M_S2_trial_varT = nanmean(nanmean(statsAll.S2.mean_trial(:,12:22,subs_var)));
    M_S1_time_conN = nanmean(nanmean(statsAll.S1.mean_time(1,:,subs_con)));
    M_S1_time_conT = nanmean(nanmean(statsAll.S1.mean_time(2,:,subs_con)));
    M_S2_time_conN = nanmean(nanmean(statsAll.S2.mean_time(1,:,subs_con)));
    M_S2_time_conT = nanmean(nanmean(statsAll.S2.mean_time(2,:,subs_con)));
    M_S1_time_varN = nanmean(nanmean(statsAll.S1.mean_time(1,:,subs_var)));
    M_S1_time_varT = nanmean(nanmean(statsAll.S1.mean_time(2,:,subs_var)));
    M_S2_time_varN = nanmean(nanmean(statsAll.S2.mean_time(1,:,subs_var)));
    M_S2_time_varT = nanmean(nanmean(statsAll.S2.mean_time(2,:,subs_var)));
    % calculate standard deviation -> 1 value per variable (std per
    % trial/bin per subject -> mean of all std values per subject -> mean
    % of all subjects std means)
    SD_S1_trial_conN = nanmean(nanmean(statsAll.S1.SD_trial(:,1:11,subs_con)));
    SD_S1_trial_conT = nanmean(nanmean(statsAll.S1.SD_trial(:,12:22,subs_con)));
    SD_S2_trial_conN = nanmean(nanmean(statsAll.S2.SD_trial(:,1:11,subs_con)));
    SD_S2_trial_conT = nanmean(nanmean(statsAll.S2.SD_trial(:,12:22,subs_con)));
    SD_S1_trial_varN = nanmean(nanmean(statsAll.S1.SD_trial(:,1:11,subs_var)));
    SD_S1_trial_varT = nanmean(nanmean(statsAll.S1.SD_trial(:,12:22,subs_var)));
    SD_S2_trial_varN = nanmean(nanmean(statsAll.S2.SD_trial(:,1:11,subs_var)));
    SD_S2_trial_varT = nanmean(nanmean(statsAll.S2.SD_trial(:,12:22,subs_var)));
    SD_S1_time_conN = nanmean(nanmean(statsAll.S1.SD_time(1,:,subs_con)));
    SD_S1_time_conT = nanmean(nanmean(statsAll.S1.SD_time(2,:,subs_con)));
    SD_S2_time_conN = nanmean(nanmean(statsAll.S2.SD_time(1,:,subs_con)));
    SD_S2_time_conT = nanmean(nanmean(statsAll.S2.SD_time(2,:,subs_con)));
    SD_S1_time_varN = nanmean(nanmean(statsAll.S1.SD_time(1,:,subs_var)));
    SD_S1_time_varT = nanmean(nanmean(statsAll.S1.SD_time(2,:,subs_var)));
    SD_S2_time_varN = nanmean(nanmean(statsAll.S2.SD_time(1,:,subs_var)));
    SD_S2_time_varT = nanmean(nanmean(statsAll.S2.SD_time(2,:,subs_var)));
    
    
    % plot SCR mean amplitude per trial
    % constant
    axes(a1); 
    plot(nanmean(statsAll.S1.mean_trial(:,1:11,subs_con),3),'LineWidth',linePlot,'Color',colConN); hold on; 
    plot(nanmean(statsAll.S1.mean_trial(:,12:22,subs_con),3),'LineWidth',linePlot,'Color',colConT);
    xlim(x_range_trial); ylim(y_range_trial); hold off;
    title({sprintf('BayBP conditioning constant') sprintf('N = %.0f', length(subs_con))},'FontWeight','bold','FontSize',12);
    xlabel(xLabel_plots_trial,'FontWeight','bold','FontSize',10);
    ylabel(yLabel_plots,'FontWeight','bold','FontSize',10); 
    legend(legend_label);
    text(x,y,sprintf('M %.3f\nSD %.3f', M_S1_trial_conN, SD_S1_trial_conN),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colConN,...
        'FontWeight','bold',...
        'FontSize',12);
    text(x+0.3,y,sprintf('M %.3f\nSD %.3f', M_S1_trial_conT, SD_S1_trial_conT),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colConT,...
        'FontWeight','bold',...
        'FontSize',12);
    box off;
    axes(a2);
    % no TENS
    B1 = bar(1,M_S1_trial_conN);
    set(B1,'FaceColor',colConN);
    hold on;
    grid on;
    errorbar(1,M_S1_trial_conN,SD_S1_trial_conN,'k','LineWidth',Line_with); % standard error mean
    % TENS
    B2 = bar(3,M_S1_trial_conT);
    set(B2,'FaceColor',colConT);
    errorbar(3,M_S1_trial_conT,SD_S1_trial_conT,'k','LineWidth',Line_with); % standard error mean
    ylim(y_range_trial);
    Labels = {'no TENS', 'TENS'};
    set(gca, 'XTick', [1 3], 'XTickLabel', Labels,'FontWeight','bold','FontSize',10);
%     [h,p,ci,stats] = ttest(MsubsStat_S1conNData,MsubsStat_S1conTData);
%     text(x+.45,y,sprintf('t(%2d)=%.2f\np=%.3f', stats.df,stats.tstat,p),...
%         'Units', 'normalized', ...
%         'VerticalAlignment', 'top', ...
%         'HorizontalAlignment', 'right', ...
%         'Color', 'k',...
%         'FontWeight','bold',...
%         'FontSize',10);
    box off;
    hold off;
    
    axes(a3);
    plot(nanmean(statsAll.S2.mean_trial(:,1:11,subs_con),3),'LineWidth',linePlot,'Color',colConN); hold on; 
    plot(nanmean(statsAll.S2.mean_trial(:,12:22,subs_con),3),'LineWidth',linePlot,'Color',colConT);
    xlim(x_range_trial); ylim(y_range_trial); hold off;
    title({sprintf('BayBP test phase constant') sprintf('N = %.0f', length(subs_con))},'FontWeight','bold','FontSize',12);
    xlabel(xLabel_plots_trial,'FontWeight','bold','FontSize',10);
    ylabel(yLabel_plots,'FontWeight','bold','FontSize',10); 
    legend(legend_label);
    text(x,y,sprintf('M %.3f\nSD %.3f', M_S2_trial_conN, SD_S2_trial_conN),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colConN,...
        'FontWeight','bold',...
        'FontSize',12);
    text(x+0.3,y,sprintf('M %.3f\nSD %.3f', M_S2_trial_conT, SD_S2_trial_conT),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colConT,...
        'FontWeight','bold',...
        'FontSize',12);
    box off;
    axes(a4);
    % no TENS
    B1 = bar(1,M_S2_trial_conN);
    set(B1,'FaceColor',colConN);
    hold on;
    grid on;
    errorbar(1,M_S2_trial_conN,SD_S2_trial_conN,'k','LineWidth',Line_with); % standard error mean
    % TENS
    B2 = bar(3,M_S2_trial_conT);
    set(B2,'FaceColor',colConT);
    errorbar(3,M_S2_trial_conT,SD_S2_trial_conT,'k','LineWidth',Line_with); % standard error mean
    ylim(y_range_trial);
    Labels = {'no TENS', 'TENS'};
    set(gca, 'XTick', [1 3], 'XTickLabel', Labels,'FontWeight','bold','FontSize',10);
    box off;
    hold off;
    % variable
    axes(a5);
    plot(nanmean(statsAll.S1.mean_trial(:,1:11,subs_var),3),'LineWidth',linePlot,'Color',colVarN); hold on; 
    plot(nanmean(statsAll.S1.mean_trial(:,12:22,subs_var),3),'LineWidth',linePlot,'Color',colVarT);
    xlim(x_range_trial); ylim(y_range_trial); hold off;
    title({sprintf('BayBP conditioning variable') sprintf('N = %.0f', length(subs_var))},'FontWeight','bold','FontSize',12);
    xlabel(xLabel_plots_trial,'FontWeight','bold','FontSize',10);
    ylabel(yLabel_plots,'FontWeight','bold','FontSize',10); 
    legend(legend_label);
    text(x,y,sprintf('M %.3f\nSD %.3f', M_S1_trial_varN, SD_S1_trial_varN),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colVarN,...
        'FontWeight','bold',...
        'FontSize',12);
    text(x+0.3,y,sprintf('M %.3f\nSD %.3f', M_S1_trial_varT, SD_S1_trial_varT),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colVarT,...
        'FontWeight','bold',...
        'FontSize',12);
    box off;
    axes(a6);
    % no TENS
    B1 = bar(1,M_S1_trial_varN);
    set(B1,'FaceColor',colVarN);
    hold on;
    grid on;
    errorbar(1,M_S1_trial_varN,SD_S1_trial_varN,'k','LineWidth',Line_with); % standard error mean
    % TENS
    B2 = bar(3,M_S1_trial_varT);
    set(B2,'FaceColor',colVarT);
    errorbar(3,M_S1_trial_varT,SD_S1_trial_varT,'k','LineWidth',Line_with); % standard error mean
    ylim(y_range_trial);
    Labels = {'no TENS', 'TENS'};
    set(gca, 'XTick', [1 3], 'XTickLabel', Labels,'FontWeight','bold','FontSize',10);
    box off;
    hold off;
    axes(a7);
    plot(nanmean(statsAll.S2.mean_trial(:,1:11,subs_var),3),'LineWidth',linePlot,'Color',colVarN); hold on; 
    plot(nanmean(statsAll.S2.mean_trial(:,12:22,subs_var),3),'LineWidth',linePlot,'Color',colVarT);
    xlim(x_range_trial); ylim(y_range_trial); hold off;
    title({sprintf('BayBP test phase variable') sprintf('N = %.0f', length(subs_var))},'FontWeight','bold','FontSize',12);
    xlabel(xLabel_plots_trial,'FontWeight','bold','FontSize',10);
    ylabel(yLabel_plots,'FontWeight','bold','FontSize',10); 
    legend(legend_label);
    text(x,y,sprintf('M %.3f\nSD %.3f', M_S2_trial_varN, SD_S2_trial_varN),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colVarN,...
        'FontWeight','bold',...
        'FontSize',12);
    text(x+0.3,y,sprintf('M %.3f\nSD %.3f', M_S2_trial_varT, SD_S2_trial_varT),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colVarT,...
        'FontWeight','bold',...
        'FontSize',12);
    box off;
    axes(a8);
    % no TENS
    B1 = bar(1,M_S2_trial_varN);
    set(B1,'FaceColor',colVarN);
    hold on;
    grid on;
    errorbar(1,M_S2_trial_varN,SD_S2_trial_varN,'k','LineWidth',Line_with); % standard error mean
    % TENS
    B2 = bar(3,M_S2_trial_varT);
    set(B2,'FaceColor',colVarT);
    errorbar(3,M_S2_trial_varT,SD_S2_trial_varT,'k','LineWidth',Line_with); % standard error mean
    ylim(y_range_trial);
    Labels = {'no TENS', 'TENS'};
    set(gca, 'XTick', [1 3], 'XTickLabel', Labels,'FontWeight','bold','FontSize',10);
    box off;
    hold off;
    
%     fname_plot = sprintf('%sBayBP_SCR_descrStatsAll_heatTrials.pdf',pathHome);
%     SaveFigure(fname_plot);
%     close gcf; 
    clear a1 a2 a3 a4 a5 a6 a7 a8 f1;
    
    % plot SCR mean amplitude per timebin
    f2 = figure;
    set(f2, 'PaperUnits','centimeters');
    set(f2, 'PaperSize',[X Y]);
    set(f2, 'PaperPosition',[xMargin yMargin xSize ySize]);
    set(f2, 'PaperOrientation','Portrait');
    b1 = axes('position', [X_left1 Y_top width_fig1 height_fig]); 
    b2 = axes('position', [X_left2 Y_top width_fig2 height_fig]);
    % subplot 2 - constant test phase
    b3 = axes('position', [X_right1 Y_top width_fig1 height_fig]);
    b4 = axes('position', [X_right2 Y_top width_fig2 height_fig]);
    % subplot 3 - variable conditioning
    b5 = axes('position', [X_left1 Y_bottom width_fig1 height_fig]);
    b6 = axes('position', [X_left2 Y_bottom width_fig2 height_fig]);
    % subplot 4 - variable test phase
    b7 = axes('position', [X_right1 Y_bottom width_fig1 height_fig]);
    b8 = axes('position', [X_right2 Y_bottom width_fig2 height_fig]);
    % constant
    axes(b1);
    plot(nanmean(statsAll.S1.mean_time(1,:,subs_con),3),'LineWidth',linePlot,'Color',colConN); hold on; 
    plot(nanmean(statsAll.S1.mean_time(2,:,subs_con),3),'LineWidth',linePlot,'Color',colConT);
    xlim(x_range_time); ylim(y_range_time); hold off;
    title({sprintf('BayBP conditioning constant') sprintf('N = %.0f', length(subs_con))},'FontWeight','bold','FontSize',12);
    xlabel(xLabel_plots_time,'FontWeight','bold','FontSize',10);
    ylabel(yLabel_plots,'FontWeight','bold','FontSize',10); 
    legend(legend_label);
    text(x,y,sprintf('M %.3f\nSD %.3f', M_S1_time_conN, SD_S1_time_conN),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colConN,...
        'FontWeight','bold',...
        'FontSize',12);
    text(x+0.3,y,sprintf('M %.3f\nSD %.3f', M_S1_time_conT, SD_S1_time_conT),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colConT,...
        'FontWeight','bold',...
        'FontSize',12);
    box off;
    axes(b2);
    % no TENS
    B1 = bar(1,M_S1_time_conN);
    set(B1,'FaceColor',colConN);
    hold on;
    grid on;
    errorbar(1,M_S1_time_conN,SD_S1_time_conN,'k','LineWidth',Line_with); % standard error mean
    % TENS
    B2 = bar(3,M_S1_time_conT);
    set(B2,'FaceColor',colConT);
    errorbar(3,M_S1_time_conT,SD_S1_time_conT,'k','LineWidth',Line_with); % standard error mean
    ylim(y_range_time);
    Labels = {'no TENS', 'TENS'};
    set(gca, 'XTick', [1 3], 'XTickLabel', Labels,'FontWeight','bold','FontSize',10);
    
    axes(b3);
    plot(nanmean(statsAll.S2.mean_time(1,:,subs_con),3),'LineWidth',linePlot,'Color',colConN); hold on; 
    plot(nanmean(statsAll.S2.mean_time(2,:,subs_con),3),'LineWidth',linePlot,'Color',colConT);
    xlim(x_range_time); ylim(y_range_time); hold off;
    title({sprintf('BayBP test phase constant') sprintf('N = %.0f', length(subs_con))},'FontWeight','bold','FontSize',12);
    xlabel(xLabel_plots_time,'FontWeight','bold','FontSize',10);
    ylabel(yLabel_plots,'FontWeight','bold','FontSize',10); 
    legend(legend_label);
    text(x,y,sprintf('M %.3f\nSD %.3f', M_S2_time_conN, SD_S2_time_conN),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colConN,...
        'FontWeight','bold',...
        'FontSize',12);
    text(x+0.3,y,sprintf('M %.3f\nSD %.3f', M_S2_time_conT, SD_S2_time_conT),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colConT,...
        'FontWeight','bold',...
        'FontSize',12);
    box off;
    axes(b4);
    % no TENS
    B1 = bar(1,M_S2_time_conN);
    set(B1,'FaceColor',colConN);
    hold on;
    grid on;
    errorbar(1,M_S2_time_conN,SD_S2_time_conN,'k','LineWidth',Line_with); % standard error mean
    % TENS
    B2 = bar(3,M_S2_time_conT);
    set(B2,'FaceColor',colConT);
    errorbar(3,M_S2_time_conT,SD_S2_time_conT,'k','LineWidth',Line_with); % standard error mean
    ylim(y_range_time);
    Labels = {'no TENS', 'TENS'};
    set(gca, 'XTick', [1 3], 'XTickLabel', Labels,'FontWeight','bold','FontSize',10);
    % variable
    axes(b5);
    plot(nanmean(statsAll.S1.mean_time(1,:,subs_var),3),'LineWidth',linePlot,'Color',colVarN); hold on; 
    plot(nanmean(statsAll.S1.mean_time(2,:,subs_var),3),'LineWidth',linePlot,'Color',colVarT);
    xlim(x_range_time); ylim(y_range_time); hold off;
    title({sprintf('BayBP conditioning variable') sprintf('N = %.0f', length(subs_var))},'FontWeight','bold','FontSize',12);
    xlabel(xLabel_plots_time,'FontWeight','bold','FontSize',10);
    ylabel(yLabel_plots,'FontWeight','bold','FontSize',10); 
    legend(legend_label);
    text(x,y,sprintf('M %.3f\nSD %.3f', M_S1_time_varN, SD_S1_time_varN),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colVarN,...
        'FontWeight','bold',...
        'FontSize',12);
    text(x+0.3,y,sprintf('M %.3f\nSD %.3f', M_S1_time_varT, SD_S1_time_varT),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colVarT,...
        'FontWeight','bold',...
        'FontSize',12);
    box off;
    axes(b6);
    % no TENS
    B1 = bar(1,M_S1_time_varN);
    set(B1,'FaceColor',colVarN);
    hold on;
    grid on;
    errorbar(1,M_S1_time_varN,SD_S1_time_varN,'k','LineWidth',Line_with); % standard error mean
    % TENS
    B2 = bar(3,M_S1_time_varT);
    set(B2,'FaceColor',colVarT);
    errorbar(3,M_S1_time_varT,SD_S1_time_varT,'k','LineWidth',Line_with); % standard error mean
    ylim(y_range_time);
    Labels = {'no TENS', 'TENS'};
    set(gca, 'XTick', [1 3], 'XTickLabel', Labels,'FontWeight','bold','FontSize',10);
    axes(b7);
    plot(nanmean(statsAll.S2.mean_time(1,:,subs_var),3),'LineWidth',linePlot,'Color',colVarN); hold on; 
    plot(nanmean(statsAll.S2.mean_time(2,:,subs_var),3),'LineWidth',linePlot,'Color',colVarT);
    xlim(x_range_time); ylim(y_range_time); hold off;
    title({sprintf('BayBP test phase variable') sprintf('N = %.0f', length(subs_var))},'FontWeight','bold','FontSize',12);
    xlabel(xLabel_plots_time,'FontWeight','bold','FontSize',10);
    ylabel(yLabel_plots,'FontWeight','bold','FontSize',10);
    legend(legend_label);
    text(x,y,sprintf('M %.3f\nSD %.3f', M_S2_time_varN, SD_S2_time_varN),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colVarN,...
        'FontWeight','bold',...
        'FontSize',12);
    text(x+0.3,y,sprintf('M %.3f\nSD %.3f', M_S2_time_varT, SD_S2_time_varT),...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', ...
        'Color', colVarT,...
        'FontWeight','bold',...
        'FontSize',12);
    box off;
    axes(b8);
    % no TENS
    B1 = bar(1,M_S2_time_varN);
    set(B1,'FaceColor',colVarN);
    hold on;
    grid on;
    errorbar(1,M_S2_time_varN,SD_S2_time_varN,'k','LineWidth',Line_with); % standard error mean
    % TENS
    B2 = bar(3,M_S2_time_varT);
    set(B2,'FaceColor',colVarT);
    errorbar(3,M_S2_time_varT,SD_S2_time_varT,'k','LineWidth',Line_with); % standard error mean
    ylim(y_range_time);
    Labels = {'no TENS', 'TENS'};
    set(gca, 'XTick', [1 3], 'XTickLabel', Labels,'FontWeight','bold','FontSize',10);
    
%     fname_plot = sprintf('%sBayBP_SCR_descrStatsAll_heatTime.pdf',pathHome);
%     SaveFigure(fname_plot);
%     close gcf;




%% plot descriptive stats
elseif strcmp(whattodo,'plot_stats')
    
    session         = varargin{1};
    subject         = varargin{2};
    target_timebin  = varargin{3};
    time_window_sec = varargin{4};
    for subject = varargin{2}
        clear descr_stats;
        load(sprintf('%stoLedalab_BayBP_%02d_descrStats.mat',path,subject));
        for session = varargin{1}
            fprintf('Processing subject #%03d... session #%d...\n',subject,session);
            close all;
            plotSizeVer = 4; % subplot size vertical
            plotSizeHor = 2; % subplot size horizontal
            X           = 25;               % x-axis paper size
            Y           = 34;               % y-axis paper size
            xMargin     = .1;               % left/right margins from page borders
            yMargin     = .1;               % bottom/top margins from page borders
            xSize       = X - 2*xMargin;    % figure size on paper (width & height)
            ySize       = Y - 2*yMargin;    % figure size on paper (width & height)
            
            F = figure;
            
            set(F, 'PaperUnits','centimeters');
            set(F, 'PaperSize',[X Y]);
            set(F, 'PaperPosition',[xMargin yMargin xSize ySize]);
            
            % session 1 - conditioning phase
            subplot(plotSizeVer,plotSizeHor,1);
            plot(descr_stats.S1.mean_time(:,1),'r'); hold on; plot(descr_stats.S1.mean_time(:,2));
            xlim([0 181]);
            legend('no TENS','TENS','Orientation','horizontal'); legend('boxoff');
            xlabel(sprintf('time (in %03d ms)',(target_timebin/time_window_sec)*10),'FontWeight','bold');
            ylabel('SCR mean amplitude (in \muS)','FontWeight','bold');
            text(.3,1.4,sprintf('BayBP%02d: SCR phasic driver - heat onset (epoch: %02d sec.)',subject, time_window_sec), ...
                'Units', 'normalized', ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'left', ...
                'Color', 'black',...
                'FontWeight','bold',...
                'FontSize',15);
            text(.8,1.2,'Session 1 - conditioning phase', ...
                'Units', 'normalized', ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'left', ...
                'Color', 'black',...
                'FontWeight','bold',...
                'FontSize',12);
            subplot(plotSizeVer,plotSizeHor,2);
            plot(descr_stats.S1.mean_trial(1:11),'r'); hold on; plot(descr_stats.S1.mean_trial(12:22));
            xlim([0 12]);
            legend('no TENS','TENS','Orientation','horizontal'); legend('boxoff');
            xlabel('trial','FontWeight','bold');
            
            subplot(plotSizeVer,plotSizeHor,3);
            plot(descr_stats.S1.variance_time(:,1),'r'); hold on; plot(descr_stats.S1.variance_time(:,2));
            xlim([0 181]);
            legend('no TENS','TENS','Orientation','horizontal'); legend('boxoff');
            xlabel(sprintf('time (in %03d ms)',(target_timebin/time_window_sec)*10),'FontWeight','bold');
            ylabel('SCR variance','FontWeight','bold');
            subplot(plotSizeVer,plotSizeHor,4);
            plot(descr_stats.S1.variance_trial(1:11),'r'); hold on; plot(descr_stats.S1.variance_trial(12:22));
            xlim([0 12]);
            legend('no TENS','TENS','Orientation','horizontal'); legend('boxoff');
            xlabel('trial','FontWeight','bold');
            
            % session 2 - placebo test phase
            subplot(plotSizeVer,plotSizeHor,5);
            plot(descr_stats.S2.mean_time(:,1),'r'); hold on; plot(descr_stats.S2.mean_time(:,2));
            xlim([0 181]);
            legend('no TENS','TENS','Orientation','horizontal'); legend('boxoff');
            xlabel(sprintf('time (in %03d ms)',(target_timebin/time_window_sec)*10),'FontWeight','bold');
            ylabel('SCR mean amplitude (in \muS)','FontWeight','bold');
            text(.9,1.2,'Session 2 - test phase', ...
                'Units', 'normalized', ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'left', ...
                'Color', 'black',...
                'FontWeight','bold',...
                'FontSize',12);
            subplot(plotSizeVer,plotSizeHor,6);
            plot(descr_stats.S2.mean_trial(1:11),'r'); hold on; plot(descr_stats.S2.mean_trial(12:22));
            xlim([0 12]);
            legend('no TENS','TENS','Orientation','horizontal'); legend('boxoff');
            xlabel('trial','FontWeight','bold');
            
            subplot(plotSizeVer,plotSizeHor,7);
            plot(descr_stats.S2.variance_time(:,1),'r'); hold on; plot(descr_stats.S2.variance_time(:,2));
            xlim([0 181]);
            legend('no TENS','TENS','Orientation','horizontal'); legend('boxoff');
            xlabel(sprintf('time (in %03d ms)',(target_timebin/time_window_sec)*10),'FontWeight','bold');
            ylabel('SCR variance','FontWeight','bold');
            subplot(plotSizeVer,plotSizeHor,8);
            plot(descr_stats.S2.variance_trial(1:11),'r'); hold on; plot(descr_stats.S2.variance_trial(12:22));
            xlim([0 12]);
            legend('no TENS','TENS','Orientation','horizontal'); legend('boxoff');
            xlabel('trial','FontWeight','bold');
        end
        
    fname_plot = sprintf('%stoLedalab_BayBP_%02d_descrStats_plot.pdf',path,subject);
    SaveFigure(fname_plot);   
        
    end
    
end


    


% bay_scr2ledalab('run',1:2,[1:3 5:19]);
% clear all; close all;
% bay_scr2ledalab('load_results',1:2,[1:3 5:19]);
% clear all; close all;
% bay_scr2ledalab('plot',1:2,[1:3 5:19]);
% clear all; close all;
% bay_scr2ledalab('get_timecourse',1:2,[1:3 5:19]);
% clear all; close all;


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
