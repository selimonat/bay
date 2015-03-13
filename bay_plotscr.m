function bay_plotscr(session,subject)


path = '/Users/onat/Desktop/toLedalab/';
load(sprintf('%stoLedalab_BayBP_S%d_%02d.mat',path,session,subject))

figure(1)
clf
plot(data.time,data.conductance)
hold on
et = [data.event(:).time];
plot([et;et],ylim,'k')
hold off;