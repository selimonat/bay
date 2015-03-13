
%%
clear all
for session = 1:2;
    for subject = 2:19
        path         =  '/Users/onat/Desktop/toLedalab/';
        fname        = sprintf('%stoLedalab_BayBP_S%d_%02d.mat',path,session,subject);
        fname_save   = sprintf('%s%s/toLedalab_BayBP_S%d_%02d.mat',path,'corrected',session,subject);
        if exist(fname_save) == 0;
            mkdir(fileparts(fname_save));
        end
        load(fname);
        for e = 1:length(data.event)
            nid  = data.event(e).nid;
            name = data.event(e).name;
            if nid == 24%cue+
                data.event(e).nid = 3;%tens/cue/stim
            elseif nid == 44%heat+
                data.event(e).nid = 4;
            elseif nid == 23%cue-
                data.event(e).nid = 1;
            elseif nid == 43%heat-
                data.event(e).nid = 2;
            end
            if strcmp(name,'W')
                data.event(e).name = '+/heat';
            elseif strcmp(name,'T')
                data.event(e).name = '+/cue'
            elseif strcmp(name,'H')
                data.event(e).name = '-/heat'
            elseif strcmp(name,'N')
                data.event(e).name = '-/cue';
            end
        end
        save(fname_save,'data');
    end
end
%%
path         =  '/Users/onat/Dropbox/Bay_SCR_Ladalab/';
for session      = 1;
    for subject      = 3
            try
                fname        = sprintf('%stoLedalab_BayBP_S%d_%02d.mat',path,session,subject);
                Ledalab({fname}, 'open', 'mat','analyze','CDA', 'optimize',6, 'overview',  1, 'export_era', [-1 15 0 1], 'export_scrlist', [0 1], 'export_eta', 1);
            end
    end
end
%%
clear all;
path         =  '/Users/onat/Desktop/toLedalab/corrected/';
session = 1;
for subject = 2:19;
    try
        fname   = sprintf('%stoLedalab_BayBP_S%d_%02d_results.mat',path,session,subject);
        load(fname);
        d = analysis.split_driver
        figure;
        plot(d.x(:,1:4),d.mean)
    end
end
%%

%%
%color code for tens vs. no-tens
co{[43]} = 'r';
co{[44]} = 'b';
co{[23]} = 'r';
co{[24]} = 'b';
%
t{[1]}   = { 'cue'};
t{[21]}  = { 'heat'};
c        = 0;
%%
figure(100);
c    = 0;
for cue = [ 0 20];
    c = c + 1;
    subplot(1,2,c);
    for block = [23 24]+cue
        i = analysis.split_phasicData.c == block ;
        plot(mean(analysis.split_phasicData.y(:,i),2),'color',co{block})
        hold on;
    end
    [cue block]
    title(t{cue+1});
    legend('tens-','tens+');
    hold off
end