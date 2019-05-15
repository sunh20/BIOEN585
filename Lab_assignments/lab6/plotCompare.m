function [] = plotCompare(est,flag)
global dat
% flag = 1 if comparing
% flag = 0 if just plotting data

if ~exist('flag','var')
    flag = 1;
end

figure;
h1 = subplot(1,4,1);
if flag == 1
    barh([dat.f300,est(:,1)])
else
    barh(dat.f300)
end
set(h1, 'Ydir', 'reverse')
yticklabels(dat.labels)
title('300 pN/s loading rate')
ylabel('Force rupture bin')
xlabel('Occurences')

h1 = subplot(1,4,2);
if flag == 1
    barh([dat.f3000,est(:,2)])
else
    barh(dat.f3000)
end
set(h1, 'Ydir', 'reverse')
yticklabels(dat.labels)
title('3000 pN/s loading rate')
xlim([0 30])

h1 = subplot(1,4,3);
if flag == 1
    barh([dat.f30000,est(:,3)])
else
    barh(dat.f30000)
end
set(h1, 'Ydir', 'reverse')
yticklabels(dat.labels)
title('30000 pN/s loading rate')
xlim([0 30])

h1 = subplot(1,4,4);
if flag == 1
    barh([dat.ctrl,est(:,4)])
else
    barh(dat.ctrl)
end
set(h1, 'Ydir', 'reverse')
yticklabels(dat.labels)
title('Negative Control')

if flag == 1
    legend('Data','Estimate')
end

if flag == 1
    figure;
    subplot(4,1,1)
    plot(dat.bins,est(:,1))
    hold on
    plot(dat.bins,dat.f300,'o')
    xlabel('rupture force')
    ylabel('occurences')
    title('300 pN/s loading rate')
    
    subplot(4,1,2)
    plot(dat.bins,est(:,2))
    hold on
    plot(dat.bins,dat.f3000,'o')
    xlabel('rupture force')
    ylabel('occurences')
    title('3000 pN/s loading rate')
    
    subplot(4,1,3)
    plot(dat.bins,est(:,3))
    hold on
    plot(dat.bins,dat.f30000,'o')
    xlabel('rupture force')
    ylabel('occurences')
    title('30000 pN/s loading rate')
    
    subplot(4,1,4)
    plot(dat.bins,est(:,4))
    hold on
    plot(dat.bins,dat.ctrl,'o')
    xlabel('rupture force')
    ylabel('occurences')
    title('Negative Control')
    
end

end