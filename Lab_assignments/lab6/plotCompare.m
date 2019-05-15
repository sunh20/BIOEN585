function [] = plotCompare(est,flag)
% plot 1: bar graph
% plot 2: line graph comparison
% plot 3: residual line graph

global dat
% flag = 1 if comparing
% flag = 0 if just plotting data

if ~exist('flag','var')
    flag = 1;
end

ref = [dat.f300, dat.f3000, dat.f30000, dat.ctrl];
name = ["300 pN/s","3000 pN/s","30000 pN/s","Neg control"];

figure;

for i = 1:size(ref,2)
    h1 = subplot(1,4,i);
    if flag == 1
        barh([ref(:,i),est(:,i)])
    else
        barh(ref(:,i))
    end
    set(h1, 'Ydir', 'reverse')
    yticklabels(dat.labels)
    title(name(i))
    if i == 1
        ylabel('Force rupture bin')
    end
    xlabel('Occurences')
end

if flag == 1
    legend('Data','Estimate')
end

% line plot comparison
if flag == 1
    figure;
    for i = 1:size(ref,2)
        % line comparison
        subplot(4,1,i);
        plot(dat.bins,est(:,i),dat.bins,ref(:,i),'o')
        if i == 4
            xlabel('rupture force')
        end
        ylabel('occurences')     
        title(name(i))
    end
end

% plot residuals
figure;
if flag == 1
    plot(dat.bins,abs(ref-est))
    xlabel('force rupture')
    ylabel('absolute residual')
    title('absolute residuals per condition')
    legend(name)
end
end