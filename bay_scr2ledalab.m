function out = bay_scr2ledalab(whattodo,varargin)

% whattodo: 
%   'run' - runs Ledalab analysis for chosen sessions & subjects CAVE:
%           define Ledalab settings before
%   'load_results' - loads results data
%   'plot' - produces .png of results data (phasic, tonic, driver, impulse 
%            response funct., condition information)
%   'get_timecourse' - produces binning matrix per subject 3 dimensional
%                      (defined timepoints (target_timebin), trials, 
%                      amplitude)

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
        SaveFigure(fname)
    end
    
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
            M       = out.analysis.split_driver.y;
            c       = out.analysis.split_driver.c;
            M       = M(:,ismember(c,[2 4]));
            %create a binning matrix
            [BM,centers]=BinningMatrix(size(M,1),target_timebin);
            Mbinned(:,:,counter)       = BM'*M;
        end
    end
    fname = sprintf('%stoLedalab_BayBP_S%d_%02d_matrix.mat',path,session,subject);
    out = Mbinned;
    save(fname, 'out');

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
