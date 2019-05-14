%This script loads all the data for the 6-OHDA headphones experiments and
%produces plots for pitch shift with bootstrapping for error bars. Written
%by Varun Saravanan, October, 2018.

clear
close all
clc

% neg_shift.birds = {'bl36yw96';'bk9gy9';'lb160rd190'};
% pos_shift.birds = {'bl47yw80';'lb199gy78';'gy46pu6';'bl26bl27'};
% no_shift.birds = {'bl20gr152';'pu63pu64';'bl23lb23';'pu81pu82';'bk21bk22';'lb140yw4';'bk24gy24';'lb24or124'};

%Can specify which birds in the set of all bird names are to be included in generating figures.
lesioned_birds = 1; %Toggle between 6-OHDA lesioned birds (set to 1) and unlesioned birds (set to 0).

if lesioned_birds
    birds_to_use = 25:40; %For 6-OHDA lesioned birds
    neg_shift_birds = [1:3 16];
    pos_shift_birds = 4:7;
    no_shift_birds = 8:15;
    %Set these:
    num_base_days = 3;
    num_shift_days = 14;
    num_days = num_base_days + num_shift_days;
    max_syls = 12;
else
    birds_to_use = 41:46; %For unlesioned birds
    neg_shift_birds = [1 4 6];
    pos_shift_birds = [2 3 5];
    no_shift_birds = [];
    %Set these:
    num_base_days = 0;
    num_shift_days = 14;
    num_days = num_base_days + num_shift_days;
    max_syls = 9;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    P(counter) = load_bird_params_2018(i);
    counter = counter + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating mean pitch change over groups - targeted syllables, same type
%syllables and different syllables. This uses the mean pitch shift of each
%syllable and averages over all syllables.

mean_neg_pitch_shift = zeros(num_base_days + num_shift_days,1);
mean_pos_pitch_shift = zeros(num_base_days + num_shift_days,1);
mean_no_pitch_shift = zeros(num_base_days + num_shift_days,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    temp2 = [];
    for j = neg_shift_birds
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_neg_pitch_shift(i) = mean(temp2);
    
    temp2 = [];
    for j = pos_shift_birds
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_pos_pitch_shift(i) = mean(temp2);
    
    temp2 = [];
    for j = no_shift_birds
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_no_pitch_shift(i) = mean(temp2);
    
end

%Sanity check: (Looks right! Leave commented)
% figure (1)
% hold all
% plot(mean_neg_pitch_shift,'r')
% plot(mean_pos_pitch_shift,'b')
% plot(mean_no_pitch_shift,'k')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To calculate the error bars for the means (SEMs), we use bootstrapping.
%Note that this implemenation differs substantially from our previous
%version. Here, we first produce a matrix of bootstrapped values from which
%we resample. Each sample gets its own Hz value of pitch from which the
%semitone values are generated. Therefore, this way explicitly includes the
%error in estimation of the mean during baseline. 1 SD of the bootstrapped
%sample gives an accurate estimate of the SEM of the dataset.

%Use these for debugging with Gordon:
clear
close all
load('workspace_variables_for_gordon_n_7.mat') %Make sure to set the right path

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.

bootstrapping_matrix = zeros(nboot,num_days,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(P(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days      %Over all days
            if lesioned_birds
                temp = P(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
            else
                temp = P(i).data{k,j};
            end
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
        if lesioned_birds
            for n = 1:nboot
                temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
                bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            end
        end
    end
end

%Now do the actual bootstrapping: 
nboot2 = 300; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,num_days);
bootstats2 = zeros(nboot2,num_days);
bootstats3 = zeros(nboot2,num_days);
for n = 1:nboot2
    
    temp_neg_birds = datasample(neg_shift_birds,length(neg_shift_birds));
    temp_data = [];
    for t = 1:length(temp_neg_birds)
        temp_syls = datasample(1:size(P(temp_neg_birds(t)).data,2),size(P(temp_neg_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_neg_birds(t)));
        end
    end
    bootstats1(n,:) = nanmean(temp_data,1);
    
    temp_pos_birds = datasample(pos_shift_birds,length(pos_shift_birds));
    temp_data = [];
    for t = 1:length(temp_pos_birds)
        temp_syls = datasample(1:size(P(temp_pos_birds(t)).data,2),size(P(temp_pos_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_pos_birds(t)));
        end
    end
    bootstats2(n,:) = nanmean(temp_data,1);
    
    if lesioned_birds
        temp_no_birds = datasample(no_shift_birds,length(no_shift_birds));
        temp_data = [];
        for t = 1:length(temp_no_birds)
            temp_syls = datasample(1:size(P(temp_no_birds(t)).data,2),size(P(temp_no_birds(t)).data,2));
            for s = 1:length(temp_syls)
                temp_pulls = datasample(1:nboot,nboot);
                temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_no_birds(t)));
            end
        end
        bootstats3(n,:) = nanmean(temp_data,1);
    end
end

std_neg_pitch_shift = nanstd(bootstats1,0,1);
std_pos_pitch_shift = nanstd(bootstats2,0,1);
if lesioned_birds
    std_no_pitch_shift = nanstd(bootstats3,0,1);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Since we are not interested in plotting baseline trends, we will zero out
%the baseline days and only plot the shift days. Additionally, we will plot
%the negative and positive shift birds in one plot and the no shift birds
%in a separate plot.

if lesioned_birds
    mean_neg_pitch_shift_no_base = mean_neg_pitch_shift(3:end);
    mean_neg_pitch_shift_no_base(1) = 0;
    mean_pos_pitch_shift_no_base = mean_pos_pitch_shift(3:end);
    mean_pos_pitch_shift_no_base(1) = 0;
    mean_no_pitch_shift_no_base = mean_no_pitch_shift(3:end);
    mean_no_pitch_shift_no_base(1) = 0;
    std_neg_pitch_shift_no_base = std_neg_pitch_shift(3:end);
    std_neg_pitch_shift_no_base(1) = 0;
    std_pos_pitch_shift_no_base = std_pos_pitch_shift(3:end);
    std_pos_pitch_shift_no_base(1) = 0;
    std_no_pitch_shift_no_base = std_no_pitch_shift(3:end);
    std_no_pitch_shift_no_base(1) = 0;
else
    mean_neg_pitch_shift_no_base = [0; mean_neg_pitch_shift];
    mean_pos_pitch_shift_no_base = [0; mean_pos_pitch_shift];
    std_neg_pitch_shift_no_base = [0 std_neg_pitch_shift];
    std_pos_pitch_shift_no_base = [0 std_pos_pitch_shift];
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 plots the change in pitch of lesioned shift birds both positive
%and negative shifts without baseline changes.

figure (1)
hold all
errorbar(0:14,mean_neg_pitch_shift_no_base,std_neg_pitch_shift_no_base,'r','lineWidth',2,'DisplayName','-1 shift')
errorbar((0:14)+0.1,mean_pos_pitch_shift_no_base,std_pos_pitch_shift_no_base,'b','lineWidth',2,'DisplayName','+1 shift')
scatter(0:14,mean_neg_pitch_shift_no_base,30,'k','fill')
scatter((0:14)+0.1,mean_pos_pitch_shift_no_base,30,'k','fill')
legend()


% 0 cents line
plot(xlim, [0 0], 'k','LineWidth',2);

ylabel('Pitch Shift from Baseline in semitones')
xlabel('Days of shift')
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2 plots drift in pitch for the equivalent time frame for no shift 
%birds.

if lesioned_birds
    figure(2)
    hold all
    errorbar(0:14,mean_no_pitch_shift_no_base,std_no_pitch_shift_no_base,'k','lineWidth',2)
    scatter(0:14,mean_no_pitch_shift_no_base,30,'k','fill')

    plot(xlim, [0 0], 'k','LineWidth',2);

    ylabel('Pitch Shift from Baseline in semitones')
    xlabel('Days of shift')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Now we need a section to do washout analysis and produce figures. The
%following section of code will do just that.

%Can specify which birds in the set of all bird names are to be included in generating figures.
% For 6-OHDA lesioned birds washout: birds_to_use = [26 40 30 31];
% neg_shift_birds = 1:2;
% pos_shift_birds = 3:4;

% For unlesioned birds washout: birds_to_use = 41:46;
% neg_shift_birds = [1 4 6];
% pos_shift_birds = [2 3 5];

clear
close all

%Set these:
lesioned_birds = 0;
washout = 'washout_full';

if strcmp(washout,'washout')
    num_base_days = 1;
else
    if lesioned_birds
        num_base_days = 17;
    else
        num_base_days = 14;
    end
end

if lesioned_birds
    birds_to_use = [26 40 30 31];
    neg_shift_birds = 1:2;
    pos_shift_birds = 3:4;
    num_shift_days = 7;
    end_of_shift = 17;
    max_syls = 12;
else
    birds_to_use = 41:46;
    neg_shift_birds = [1 4 6];
    pos_shift_birds = [2 3 5];
    num_shift_days = 6;
    end_of_shift = 1;
    max_syls = 9;
end

num_days = num_base_days + num_shift_days;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    %Use 'washout' to load only washout data. Use 'washout_full' to load
    %pitch shift and washout data.
    Q(counter) = load_bird_params_2018(i,washout); %Loading only washout data
    counter = counter + 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating mean pitch change over groups - targeted syllables, same type
%syllables and different syllables. This uses the mean pitch shift of each
%syllable and averages over all syllables.

mean_neg_pitch_shift = zeros(num_base_days + num_shift_days,1);
mean_pos_pitch_shift = zeros(num_base_days + num_shift_days,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    temp2 = [];
    for j = neg_shift_birds
        if lesioned_birds
            temp = horzcat(Q(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(Q(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_neg_pitch_shift(i) = mean(temp2);
    
    temp2 = [];
    for j = pos_shift_birds
        if lesioned_birds
            temp = horzcat(Q(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(Q(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_pos_pitch_shift(i) = mean(temp2);
    
end

%Sanity check: (Looks right! Leave commented)
% figure (1)
% hold all
% plot(mean_neg_pitch_shift,'r')
% plot(mean_pos_pitch_shift,'b')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To calculate the error bars for the means (SEMs), we use bootstrapping.
%Note that this implemenation differs substantially from our previous
%version. Here, we first produce a matrix of bootstrapped values from which
%we resample. Each sample gets its own Hz value of pitch from which the
%semitone values are generated. Therefore, this way explicitly includes the
%error in estimation of the mean during baseline. 1 SD of the bootstrapped
%sample gives an accurate estimate of the SEM of the dataset.

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.
if lesioned_birds
    num_days_bootstrap = size(Q(1).data_hz,1);
    max_syls = 12;
else
    num_days_bootstrap = size(Q(1).data,1);
    max_syls = 9;
end

bootstrapping_matrix = zeros(nboot,num_days_bootstrap,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(Q(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days_bootstrap      %Over all days
            if lesioned_birds
                temp = Q(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
            else
                if strcmp(washout,'washout')
                    temp = Q(i).data_unsubtracted{k,j};
                else
                    temp = Q(i).data{k,j};
                end
            end
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
        for n = 1:nboot
            if lesioned_birds
                temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
                bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            end
            %Subtract last day of shift from the last day onwards:
            if strcmp(washout,'washout')
                bootstrapping_matrix(n,end_of_shift:end,j,i) = bootstrapping_matrix(n,end_of_shift:end,j,i) - bootstrapping_matrix(n,end_of_shift,j,i);
            end
        end
    end
end

if strcmp(washout,'washout')
    %Redefine matrix as last day of shift till end:
    bootstrapping_matrix = bootstrapping_matrix(:,end_of_shift:end,:,:);
end

%Now do the actual bootstrapping: 
nboot2 = 300; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,num_days);
bootstats2 = zeros(nboot2,num_days);
for n = 1:nboot2
    
    temp_neg_birds = datasample(neg_shift_birds,length(neg_shift_birds));
    temp_data = [];
    for t = 1:length(temp_neg_birds)
        temp_syls = datasample(1:size(Q(temp_neg_birds(t)).data,2),size(Q(temp_neg_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_neg_birds(t)));
        end
    end
    bootstats1(n,:) = nanmean(temp_data,1);
    
    temp_pos_birds = datasample(pos_shift_birds,length(pos_shift_birds));
    temp_data = [];
    for t = 1:length(temp_pos_birds)
        temp_syls = datasample(1:size(Q(temp_pos_birds(t)).data,2),size(Q(temp_pos_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_pos_birds(t)));
        end
    end
    bootstats2(n,:) = nanmean(temp_data,1);
    
    
end

std_neg_pitch_shift = nanstd(bootstats1,0,1);
std_pos_pitch_shift = nanstd(bootstats2,0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 plots the non-adjusted washout mean and error for the two groups
%as is.

figure (1)
hold all
errorbar(0:(num_days-1),mean_neg_pitch_shift,std_neg_pitch_shift,'r','lineWidth',2,'DisplayName','-1 shift')
errorbar((0:(num_days-1))+0.1,mean_pos_pitch_shift,std_pos_pitch_shift,'b','lineWidth',2,'DisplayName','+1 shift')
scatter(0:(num_days-1),mean_neg_pitch_shift,30,'k','fill')
scatter((0:(num_days-1))+0.1,mean_pos_pitch_shift,30,'k','fill')
legend()


% 0 cents line
plot(xlim, [0 0], 'k','LineWidth',2);

ylabel('Washout pitch change in semitones')
xlabel('Days of washout')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2 plots the washout alone in washout-full mode without the change
%during pitch shift or baseline.

if strcmp(washout,'washout_full')
    figure (2)
    hold all
    errorbar(0:(num_shift_days),mean_neg_pitch_shift(num_base_days:num_days),std_neg_pitch_shift(num_base_days:num_days),'r','lineWidth',2,'DisplayName','-1 shift')
    errorbar((0:(num_shift_days))+0.1,mean_pos_pitch_shift(num_base_days:num_days),std_pos_pitch_shift(num_base_days:num_days),'b','lineWidth',2,'DisplayName','+1 shift')
    scatter(0:(num_shift_days),mean_neg_pitch_shift(num_base_days:num_days),30,'k','fill')
    scatter((0:(num_shift_days))+0.1,mean_pos_pitch_shift(num_base_days:num_days),30,'k','fill')
    legend()


    % 0 cents line
    if lesioned_birds
        xlim([-0.2 7.2])
    else
        xlim([-0.2 6.2])
    end
    plot(xlim, [0 0], 'k','LineWidth',2);

    ylabel('Washout pitch change in semitones')
    xlabel('Days of washout')
    
    ylim([-0.8 0.6])
end

%%
%Now we need a section to produce washout figures for no shift birds.

clear
close all

%Set these:
lesioned_birds = 1;
washout = 'washout';

if strcmp(washout,'washout')
    num_base_days = 1;
else
    if lesioned_birds
        num_base_days = 17;
    else
        num_base_days = 14;
    end
end

if lesioned_birds
    birds_to_use = 35:39;
    no_shift_birds = 1:5;
    num_shift_days = 7;
    end_of_shift = 17;
    max_syls = 6;
end

num_days = num_base_days + num_shift_days;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    %Use 'washout' to load only washout data. Use 'washout_full' to load
    %pitch shift and washout data.
    Q(counter) = load_bird_params_2018(i,washout); %Loading only washout data
    counter = counter + 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating mean pitch change over groups - targeted syllables, same type
%syllables and different syllables. This uses the mean pitch shift of each
%syllable and averages over all syllables.

mean_no_pitch_shift = zeros(num_base_days + num_shift_days,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    temp2 = [];
    for j = no_shift_birds
        if lesioned_birds
            temp = horzcat(Q(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(Q(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_no_pitch_shift(i) = mean(temp2);
    
end

% Uncomment the following lines to remove baseline from "washout_full":
if strcmp(washout,'washout_full')
    mean_no_pitch_shift = mean_no_pitch_shift(3:end);
    mean_no_pitch_shift(1) = 0;
end

%Sanity check: (Looks right! Leave commented)
% figure (1)
% hold all
% plot(mean_neg_pitch_shift,'r')
% plot(mean_pos_pitch_shift,'b')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To calculate the error bars for the means (SEMs), we use bootstrapping.
%Note that this implemenation differs substantially from our previous
%version. Here, we first produce a matrix of bootstrapped values from which
%we resample. Each sample gets its own Hz value of pitch from which the
%semitone values are generated. Therefore, this way explicitly includes the
%error in estimation of the mean during baseline. 1 SD of the bootstrapped
%sample gives an accurate estimate of the SEM of the dataset.

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.
if lesioned_birds
    num_days_bootstrap = size(Q(1).data_hz,1);
    max_syls = 6;
else
    num_days_bootstrap = size(Q(1).data,1);
    max_syls = 9;
end

bootstrapping_matrix = zeros(nboot,num_days_bootstrap,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(Q(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days_bootstrap      %Over all days
            if lesioned_birds
                temp = Q(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
            else
                if strcmp(washout,'washout')
                    temp = Q(i).data_unsubtracted{k,j};
                else
                    temp = Q(i).data{k,j};
                end
            end
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
        for n = 1:nboot
            if lesioned_birds
                temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
                bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            end
            %Subtract last day of shift from the last day onwards:
            if strcmp(washout,'washout')
                bootstrapping_matrix(n,end_of_shift:end,j,i) = bootstrapping_matrix(n,end_of_shift:end,j,i) - bootstrapping_matrix(n,end_of_shift,j,i);
            end
        end
    end
end

if strcmp(washout,'washout')
    %Redefine matrix as last day of shift till end:
    bootstrapping_matrix = bootstrapping_matrix(:,end_of_shift:end,:,:);
end

%Now do the actual bootstrapping: 
nboot2 = 300; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,num_days);
for n = 1:nboot2
    
    temp_no_birds = datasample(no_shift_birds,length(no_shift_birds));
    temp_data = [];
    for t = 1:length(temp_no_birds)
        temp_syls = datasample(1:size(Q(temp_no_birds(t)).data,2),size(Q(temp_no_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_no_birds(t)));
        end
    end
    bootstats1(n,:) = nanmean(temp_data,1);
    
    
end

std_no_pitch_shift = nanstd(bootstats1,0,1);

% Uncomment the following lines to remove baseline from "washout_full":
if strcmp(washout,'washout_full')
    std_no_pitch_shift = std_no_pitch_shift(3:end);
    std_no_pitch_shift(1) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 plots the non-adjusted washout mean and error for the two groups
%as is.

figure (1)
hold all
errorbar(0:(length(mean_no_pitch_shift)-1),mean_no_pitch_shift,std_no_pitch_shift,'k','lineWidth',2)
scatter(0:(length(mean_no_pitch_shift)-1),mean_no_pitch_shift,30,'k','fill')

% 0 cents line
plot(xlim, [0 0], 'k','LineWidth',2);
if strcmp(washout,'washout_full')
    plot([14.5 14.5],ylim,'k--','LineWidth',2);
end
ylabel('Washout pitch change in semitones')
xlabel('Days of washout')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2 plots the washout alone in washout-full mode without the change
%during pitch shift or baseline.

if strcmp(washout,'washout_full')
    figure (2)
    hold all
    errorbar(0:(num_shift_days),mean_no_pitch_shift((num_base_days-2):(num_days-2)),std_no_pitch_shift((num_base_days-2):(num_days-2)),'k','lineWidth',2)
    scatter(0:(num_shift_days),mean_no_pitch_shift((num_base_days-2):(num_days-2)),30,'k','fill')
    

    % 0 cents line
    if lesioned_birds
        xlim([-0.2 7.2])
    else
        xlim([-0.2 6.2])
    end
    plot(xlim, [0 0], 'k','LineWidth',2);

    ylabel('Washout pitch change in semitones')
    xlabel('Days of washout')
    
    ylim([-0.8 0.6])
end

%%
%Now we need a section to produce washout figures for no shift birds.

clear
close all

%Set these:
lesioned_birds = 1;
washout = 'washout_full';

if strcmp(washout,'washout')
    num_base_days = 1;
else
    if lesioned_birds
        num_base_days = 17;
    else
        num_base_days = 14;
    end
end

if lesioned_birds
    birds_to_use = 35:39;
    no_shift_birds = 1:5;
    num_shift_days = 7;
    end_of_shift = 17;
    max_syls = 6;
end

num_days = num_base_days + num_shift_days;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    %Use 'washout' to load only washout data. Use 'washout_full' to load
    %pitch shift and washout data.
    Q(counter) = load_bird_params_2018(i,washout); %Loading only washout data
    counter = counter + 1;
end



num_birds = length(birds_to_use);
mean_no_pitch_shift = zeros(1 + num_shift_days,num_birds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = end_of_shift:num_days
    if lesioned_birds
        for j = no_shift_birds
            temp = [];
            if lesioned_birds
                temp = horzcat(Q(j).data{i,:});
            else
                temp = vertcat(Q(j).data{i,:});
            end
            mean_no_pitch_shift((i-end_of_shift+1),j) = mean(temp);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 plots the individual birds for the lesioned no shift headphones
%group.

figure (1)
hold all
for j = no_shift_birds
    plot(0:(num_shift_days),mean_no_pitch_shift(:,j),'Color',[0.1*(j) 0.1*(j) 0.1*(j)],'LineWidth',2)
end
xlim([-0.2 7.2])
plot(xlim, [0 0], 'k','LineWidth',2);
ylim([-0.8 0.8])

ylabel('Washout Pitch shift')
xlabel('Days of washout')



%%
%Now we need a section of code that can produce the washout figures that
%shows for individual birds the last day of shift and combined with the
%period of washout.

clear
close all
clc

%Set these:
lesioned_birds = 0;
washout = 'washout_full';
plot_error_bars = 0;

num_base_days = 1;

if lesioned_birds
    birds_to_use = [26 40 30 31];
    neg_shift_birds = 1:2;
    pos_shift_birds = 3:4;
    num_shift_days = 7;
    end_of_shift = 17;
    max_syls = 12;
else
    birds_to_use = 41:46;
    neg_shift_birds = [1 4 6];
    pos_shift_birds = [2 3 5];
    num_shift_days = 6;
    end_of_shift = 1;
    max_syls = 9;
end

num_days = num_base_days + num_shift_days;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    %Use 'washout' to load only washout data. Use 'washout_full' to load
    %pitch shift and washout data.
    Q(counter) = load_bird_params_2018(i,washout); %Loading only washout data
    counter = counter + 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating mean pitch change over groups - targeted syllables, same type
%syllables and different syllables. This uses the mean pitch shift of each
%syllable and averages over all syllables.

mean_neg_pitch_shift = zeros(num_base_days + num_shift_days,length(neg_shift_birds));
mean_pos_pitch_shift = zeros(num_base_days + num_shift_days,length(pos_shift_birds));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    counter = 1;
    for j = neg_shift_birds
        if lesioned_birds
            temp = horzcat(Q(j).data_new{i,:});
        else
            temp = vertcat(Q(j).data_new{i,:});
        end
        mean_neg_pitch_shift(i,counter) = mean(temp);
        counter = counter + 1;
    end
    
    counter = 1;
    for j = pos_shift_birds
        if lesioned_birds
            temp = horzcat(Q(j).data_new{i,:});
        else
            temp = vertcat(Q(j).data_new{i,:});
        end
        mean_pos_pitch_shift(i,counter) = mean(temp);
        counter = counter + 1;
    end
    
end

%Sanity check: (Looks right! Leave commented)
% figure (1)
% hold all
% for i = 1:length(neg_shift_birds)
%     plot(mean_neg_pitch_shift(:,i),'r')
%     plot(mean_pos_pitch_shift(:,i),'b')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To calculate the error bars for the means (SEMs), we use bootstrapping.
%Note that this implemenation differs substantially from our previous
%version. Here, we first produce a matrix of bootstrapped values from which
%we resample. Each sample gets its own Hz value of pitch from which the
%semitone values are generated. Therefore, this way explicitly includes the
%error in estimation of the mean during baseline. 1 SD of the bootstrapped
%sample gives an accurate estimate of the SEM of the dataset.

%Bootstrapping part
if plot_error_bars
nboot = 1000; %No of times to resample for bootstrapping.
if lesioned_birds
    num_days_bootstrap = size(Q(1).data_new,1);
    max_syls = 12;
else
    num_days_bootstrap = size(Q(1).data_new,1);
    max_syls = 9;
end

bootstrapping_matrix = zeros(nboot,num_days_bootstrap,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(Q(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days_bootstrap      %Over all days
            temp = Q(i).data_new{k,j};
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
    end
end

%Now do the actual bootstrapping: 
nboot2 = 300; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,num_days);
bootstats2 = zeros(nboot2,num_days);
for i = 1:length(neg_shift_birds)
    for n = 1:nboot2

        temp_neg_birds = neg_shift_birds(i);
        temp_data = [];
        for t = 1:length(temp_neg_birds)
            temp_syls = datasample(1:size(Q(temp_neg_birds(t)).data,2),size(Q(temp_neg_birds(t)).data,2));
            for s = 1:length(temp_syls)
                temp_pulls = datasample(1:nboot,nboot);
                temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_neg_birds(t)));
            end
        end
        bootstats1(n,:) = nanmean(temp_data,1);

        temp_pos_birds = pos_shift_birds(i);
        temp_data = [];
        for t = 1:length(temp_pos_birds)
            temp_syls = datasample(1:size(Q(temp_pos_birds(t)).data,2),size(Q(temp_pos_birds(t)).data,2));
            for s = 1:length(temp_syls)
                temp_pulls = datasample(1:nboot,nboot);
                temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_pos_birds(t)));
            end
        end
        bootstats2(n,:) = nanmean(temp_data,1);


    end
    std_neg_pitch_shift(i,:) = nanstd(bootstats1,0,1);
    std_pos_pitch_shift(i,:) = nanstd(bootstats2,0,1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 plots the non-adjusted washout mean and error for the two groups
%as is.

figure (1)
hold all
for i = 1:length(neg_shift_birds)
    if plot_error_bars
        errorbar(0:(num_days-1),mean_neg_pitch_shift(:,i),std_neg_pitch_shift(i,:)','r','lineWidth',2)
        errorbar((0:(num_days-1))+0.1,mean_pos_pitch_shift(:,i),std_pos_pitch_shift(i,:)','b','lineWidth',2)
        scatter(0:(num_days-1),mean_neg_pitch_shift(:,i),30,'k','fill')
        scatter((0:(num_days-1))+0.1,mean_pos_pitch_shift(:,i),30,'k','fill')
    else
        plot(0:(num_days-1),mean_neg_pitch_shift(:,i),'Color',[1 0.5 0.5],'lineWidth',2)
        plot(0:(num_days-1),mean_pos_pitch_shift(:,i),'Color',[0.5 0.5 1],'lineWidth',2)
    end
end

% 0 cents line
ylim([-0.8 0.8])
if lesioned_birds
    xlim([-0.2 7.2])
else
    xlim([-0.2 6.2])
end
plot(xlim, [0 0], 'k','LineWidth',2);

ylabel('Washout pitch change in semitones')
xlabel('Days of washout')

%%
%The following sections of code are for special cases of the above general
%structure. The first is to produce a figure with no shift groups - both
%headphones and no headphones separately.

clear
close all
clc

% neg_shift.birds = {'bl36yw96';'bk9gy9';'lb160rd190'};
% pos_shift.birds = {'bl47yw80';'lb199gy78';'gy46pu6';'bl26bl27'};
% no_shift.birds = {'bl20gr152';'pu63pu64';'bl23lb23';'pu81pu82';'bk21bk22';'lb140yw4';'bk24gy24';'lb24or124'};

%Can specify which birds in the set of all bird names are to be included in generating figures.
lesioned_birds = 1; %Toggle between 6-OHDA lesioned birds (set to 1) and unlesioned birds (set to 0).

if lesioned_birds
    birds_to_use = 32:39; %For 6-OHDA lesioned birds
    no_shift_no_headp_birds = 1:3;
    no_shift_headp_birds = 4:8;
    %Set these:
    num_base_days = 3;
    num_shift_days = 14;
    num_days = num_base_days + num_shift_days;
    max_syls = 12;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    P(counter) = load_bird_params_2018(i);
    counter = counter + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating mean pitch change over groups - targeted syllables, same type
%syllables and different syllables. This uses the mean pitch shift of each
%syllable and averages over all syllables.

mean_no_pitch_shift_no_headp = zeros(num_base_days + num_shift_days,1);
mean_no_pitch_shift_headp = zeros(num_base_days + num_shift_days,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    temp2 = [];
    for j = no_shift_no_headp_birds
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_no_pitch_shift_no_headp(i) = mean(temp2);
    
    temp2 = [];
    for j = no_shift_headp_birds
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_no_pitch_shift_headp(i) = mean(temp2);
    
end

%Bootstrapping part:
nboot = 1000; %No of times to resample for bootstrapping.

bootstrapping_matrix = zeros(nboot,num_days,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(P(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days      %Over all days
            temp = P(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
        for n = 1:nboot
            temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
            bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            
        end
    end
end

%Now do the actual bootstrapping: 
nboot2 = 300; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,num_days);
bootstats2 = zeros(nboot2,num_days);
for n = 1:nboot2
    
    temp_no_headp_birds = datasample(no_shift_no_headp_birds,length(no_shift_no_headp_birds));
    temp_data = [];
    for t = 1:length(temp_no_headp_birds)
        temp_syls = datasample(1:size(P(temp_no_headp_birds(t)).data,2),size(P(temp_no_headp_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_no_headp_birds(t)));
        end
    end
    bootstats1(n,:) = nanmean(temp_data,1);
    
    temp_pos_birds = datasample(no_shift_headp_birds,length(no_shift_headp_birds));
    temp_data = [];
    for t = 1:length(temp_pos_birds)
        temp_syls = datasample(1:size(P(temp_pos_birds(t)).data,2),size(P(temp_pos_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_pos_birds(t)));
        end
    end
    bootstats2(n,:) = nanmean(temp_data,1);
    
end

std_no_headp_pitch_shift = nanstd(bootstats1,0,1);
std_headp_pitch_shift = nanstd(bootstats2,0,1);

if lesioned_birds
    mean_neg_pitch_shift_no_base = mean_no_pitch_shift_no_headp(3:end);
    mean_neg_pitch_shift_no_base(1) = 0;
    mean_pos_pitch_shift_no_base = mean_no_pitch_shift_headp(3:end);
    mean_pos_pitch_shift_no_base(1) = 0;
    std_no_headp_pitch_shift_no_base = std_no_headp_pitch_shift(3:end);
    std_no_headp_pitch_shift_no_base(1) = 0;
    std_headp_pitch_shift_no_base = std_headp_pitch_shift(3:end);
    std_headp_pitch_shift_no_base(1) = 0;
else
    mean_neg_pitch_shift_no_base = [0; mean_neg_pitch_shift];
    mean_pos_pitch_shift_no_base = [0; mean_pos_pitch_shift];
    std_neg_pitch_shift_no_base = [0; std_neg_pitch_shift];
    std_pos_pitch_shift_no_base = [0; std_pos_pitch_shift];
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 plots the change in pitch of lesioned shift birds both positive
%and negative shifts without baseline changes.

figure (1)
hold all
errorbar(0:14,mean_neg_pitch_shift_no_base,std_no_headp_pitch_shift_no_base,'Color',[0.5 0.5 0.5],'lineWidth',2,'DisplayName','-1 shift')
errorbar((0:14)+0.1,mean_pos_pitch_shift_no_base,std_headp_pitch_shift_no_base,'k','lineWidth',2,'DisplayName','+1 shift')
scatter(0:14,mean_neg_pitch_shift_no_base,30,'k','fill')
scatter((0:14)+0.1,mean_pos_pitch_shift_no_base,30,'k','fill')
legend()


% 0 cents line
plot(xlim, [0 0], 'k','LineWidth',2);

ylabel('Pitch Shift from Baseline in semitones')
xlabel('Days of shift')

%%
%The second is to look at the full graph for gr193gr194 fully since it is
%the special case bird.

clear
close all
clc

% neg_shift.birds = {'bl36yw96';'bk9gy9';'lb160rd190'};
% pos_shift.birds = {'bl47yw80';'lb199gy78';'gy46pu6';'bl26bl27'};
% no_shift.birds = {'bl20gr152';'pu63pu64';'bl23lb23';'pu81pu82';'bk21bk22';'lb140yw4';'bk24gy24';'lb24or124'};

%Can specify which birds in the set of all bird names are to be included in generating figures.
lesioned_birds = 1; %Toggle between 6-OHDA lesioned birds (set to 1) and unlesioned birds (set to 0).

if lesioned_birds
    birds_to_use = 40; %For 6-OHDA lesioned birds
    the_bird = 1;
    %Set these:
    num_base_days = 6;
    num_shift_days = 14;
    num_days = num_base_days + num_shift_days;
    max_syls = 12;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    P(counter) = load_bird_params_2018(i);
    counter = counter + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating mean pitch change over groups - targeted syllables, same type
%syllables and different syllables. This uses the mean pitch shift of each
%syllable and averages over all syllables.

mean_pitch_shift = zeros(num_base_days + num_shift_days,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    temp2 = [];
    for j = the_bird
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_pitch_shift(i) = mean(temp2);
    
end

%Bootstrapping part:
nboot = 1000; %No of times to resample for bootstrapping.

bootstrapping_matrix = zeros(nboot,num_days,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(P(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days      %Over all days
            temp = P(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
        for n = 1:nboot
            temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
            bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            
        end
    end
end

%Now do the actual bootstrapping: 
nboot2 = 300; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,num_days);
for n = 1:nboot2
    
    temp_no_headp_birds = datasample(the_bird,length(the_bird));
    temp_data = [];
    for t = 1:length(temp_no_headp_birds)
        temp_syls = datasample(1:size(P(temp_no_headp_birds(t)).data,2),size(P(temp_no_headp_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_no_headp_birds(t)));
        end
    end
    bootstats1(n,:) = nanmean(temp_data,1);
    
end

std_no_headp_pitch_shift = nanstd(bootstats1,0,1);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 plots the change in pitch of lesioned shift birds both positive
%and negative shifts without baseline changes.

figure (1)
hold all
errorbar(-5:14,mean_pitch_shift,std_no_headp_pitch_shift,'r','lineWidth',2,'DisplayName','-1 shift')
scatter(-5:14,mean_pitch_shift,30,'k','fill')
ylim([-1 0.6])

% 0 cents line
plot(xlim, [0 0], 'k','LineWidth',2);

% Baseline separator:
plot([0.5 0.5], ylim, 'k', 'LineWidth',2);

ylabel('Pitch Shift from Baseline in semitones')
xlabel('Days of shift')

%%
%The following section preps the data to get adaptive pitch changes for
%both unlesioned and lesioned birds. I ran this section twice, once for
%lesioned and once for unlesioned to save the mean and error for each group
%which was then plotted in the next section.

clear
close all
clc

% neg_shift.birds = {'bl36yw96';'bk9gy9';'lb160rd190'};
% pos_shift.birds = {'bl47yw80';'lb199gy78';'gy46pu6';'bl26bl27'};
% no_shift.birds = {'bl20gr152';'pu63pu64';'bl23lb23';'pu81pu82';'bk21bk22';'lb140yw4';'bk24gy24';'lb24or124'};

%Can specify which birds in the set of all bird names are to be included in generating figures.
lesioned_birds = 0; %Toggle between 6-OHDA lesioned birds (set to 1) and unlesioned birds (set to 0).

if lesioned_birds
    birds_to_use = 25:40; %For 6-OHDA lesioned birds
    neg_shift_birds = [1:3 16];
    pos_shift_birds = 4:7;
    no_shift_birds = 8:15;
    %Set these:
    num_base_days = 3;
    num_shift_days = 14;
    num_days = num_base_days + num_shift_days;
    max_syls = 12;
else
    birds_to_use = 41:46; %For unlesioned birds
    neg_shift_birds = [1 4 6];
    pos_shift_birds = [2 3 5];
    no_shift_birds = [];
    %Set these:
    num_base_days = 0;
    num_shift_days = 14;
    num_days = num_base_days + num_shift_days;
    max_syls = 9;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    P(counter) = load_bird_params_2018(i);
    counter = counter + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating mean pitch change over groups - targeted syllables, same type
%syllables and different syllables. This uses the mean pitch shift of each
%syllable and averages over all syllables.

mean_adap_pitch_shift = zeros(num_base_days + num_shift_days,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    temp2 = [];
    for j = neg_shift_birds
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    
    for j = pos_shift_birds
        if lesioned_birds
            temp = -1*horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = -1*vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_adap_pitch_shift(i) = mean(temp2);
    
end

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.

bootstrapping_matrix = zeros(nboot,num_days,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(P(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days      %Over all days
            if lesioned_birds
                temp = P(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
            else
                temp = P(i).data{k,j};
            end
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
        for n = 1:nboot
            if lesioned_birds
                temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
                bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            end
            bootstrapping_matrix(n,:,j,i) = -1*P(i).shift_direction*bootstrapping_matrix(n,:,j,i);
            
        end
    end
end

%Now do the actual bootstrapping: 
nboot2 = 300; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,num_days);
error_birds = [neg_shift_birds pos_shift_birds];
for n = 1:nboot2
    
    temp_neg_birds = datasample(error_birds,length(error_birds));
    temp_data = [];
    for t = 1:length(temp_neg_birds)
        temp_syls = datasample(1:size(P(temp_neg_birds(t)).data,2),size(P(temp_neg_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_neg_birds(t)));
        end
    end
    bootstats1(n,:) = nanmean(temp_data,1);
    
    
end

std_adap_pitch_shift = nanstd(bootstats1,0,1);

if lesioned_birds
    mean_adap_pitch_shift_no_base = mean_adap_pitch_shift(3:end);
    mean_adap_pitch_shift_no_base(1) = 0;
    std_adap_pitch_shift_no_base = std_adap_pitch_shift(3:end);
    std_adap_pitch_shift_no_base(1) = 0;
else
    mean_adap_pitch_shift_no_base = [0; mean_adap_pitch_shift];
    std_adap_pitch_shift_no_base = [0 std_adap_pitch_shift];
end

% figure (1)
% hold all
% errorbar(0:14,mean_adap_pitch_shift_no_base,std_adap_pitch_shift_no_base,'k','lineWidth',2,'DisplayName','-1 shift')
% scatter(0:14,mean_adap_pitch_shift_no_base,30,'k','fill')

%%
%In this section, we use the 2 files obtained above to plot the adaptive
%change for the 2 conditions in the same plot.

clear
close all
clc

lesioned = load('adap_change_lesioned_birds.mat');
unlesioned = load('adap_change_unlesioned_birds.mat');

figure (1)
hold all
errorbar(0:14,lesioned.mean_adap_pitch_shift_no_base,lesioned.std_adap_pitch_shift_no_base,'Color',[0.5 0.5 0.5],'lineWidth',2,'DisplayName','-1 shift')
errorbar(0:14,unlesioned.mean_adap_pitch_shift_no_base,unlesioned.std_adap_pitch_shift_no_base,'k','lineWidth',2,'DisplayName','-1 shift')
scatter(0:14,lesioned.mean_adap_pitch_shift_no_base,30,'k','fill')
scatter(0:14,unlesioned.mean_adap_pitch_shift_no_base,30,'k','fill')

% 0 cents line
xlim([-0.2 14.4])
plot(xlim, [0 0], 'k','LineWidth',2);

ylabel('Pitch change in semitones')
xlabel('Days of shift')

ylim([-1 0.6])

%%
%The following section is used to produce a plot that compares shift versus
%no shift for lesioned birds directly across all days.

clear
close all
clc

% neg_shift.birds = {'bl36yw96';'bk9gy9';'lb160rd190'};
% pos_shift.birds = {'bl47yw80';'lb199gy78';'gy46pu6';'bl26bl27'};
% no_shift.birds = {'bl20gr152';'pu63pu64';'bl23lb23';'pu81pu82';'bk21bk22';'lb140yw4';'bk24gy24';'lb24or124'};

%Can specify which birds in the set of all bird names are to be included in generating figures.
lesioned_birds = 1; %Toggle between 6-OHDA lesioned birds (set to 1) and unlesioned birds (set to 0).

if lesioned_birds
    birds_to_use = 25:40; %For 6-OHDA lesioned birds
    neg_shift_birds = [1:3 16];
    pos_shift_birds = 4:7;
    no_shift_birds = 8:15;
    %Set these:
    num_base_days = 3;
    num_shift_days = 14;
    num_days = num_base_days + num_shift_days;
    max_syls = 12;
else
    birds_to_use = 41:46; %For unlesioned birds
    neg_shift_birds = [1 4 6];
    pos_shift_birds = [2 3 5];
    no_shift_birds = [];
    %Set these:
    num_base_days = 0;
    num_shift_days = 14;
    num_days = num_base_days + num_shift_days;
    max_syls = 9;
end

shift_birds = [neg_shift_birds pos_shift_birds];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    P(counter) = load_bird_params_2018(i);
    counter = counter + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating mean pitch change over groups - targeted syllables, same type
%syllables and different syllables. This uses the mean pitch shift of each
%syllable and averages over all syllables.

mean_shift_pitch_shift = zeros(num_base_days + num_shift_days,1);
mean_no_pitch_shift = zeros(num_base_days + num_shift_days,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    temp2 = [];
    for j = shift_birds
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_shift_pitch_shift(i) = mean(temp2);
    
    temp2 = [];
    for j = no_shift_birds
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,:});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_no_pitch_shift(i) = mean(temp2);
    
end

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.

bootstrapping_matrix = zeros(nboot,num_days,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(P(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days      %Over all days
            if lesioned_birds
                temp = P(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
            else
                temp = P(i).data{k,j};
            end
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
        if lesioned_birds
            for n = 1:nboot
                temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
                bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            end
        end
    end
end

%Now do the actual bootstrapping: 
nboot2 = 300; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,num_days);
bootstats2 = zeros(nboot2,num_days);
for n = 1:nboot2
    
    temp_neg_birds = datasample(shift_birds,length(shift_birds));
    temp_data = [];
    for t = 1:length(temp_neg_birds)
        temp_syls = datasample(1:size(P(temp_neg_birds(t)).data,2),size(P(temp_neg_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_neg_birds(t)));
        end
    end
    bootstats1(n,:) = nanmean(temp_data,1);
    
    if lesioned_birds
        temp_no_birds = datasample(no_shift_birds,length(no_shift_birds));
        temp_data = [];
        for t = 1:length(temp_no_birds)
            temp_syls = datasample(1:size(P(temp_no_birds(t)).data,2),size(P(temp_no_birds(t)).data,2));
            for s = 1:length(temp_syls)
                temp_pulls = datasample(1:nboot,nboot);
                temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,:,temp_syls(s),temp_no_birds(t)));
            end
        end
        bootstats2(n,:) = nanmean(temp_data,1);
    end
end

std_pitch_shift = nanstd(bootstats1,0,1);
if lesioned_birds
    std_no_pitch_shift = nanstd(bootstats2,0,1);
end

if lesioned_birds
    mean_shift_pitch_shift_no_base = mean_shift_pitch_shift(3:end);
    mean_shift_pitch_shift_no_base(1) = 0;
    mean_no_pitch_shift_no_base = mean_no_pitch_shift(3:end);
    mean_no_pitch_shift_no_base(1) = 0;
    std_pitch_shift_no_base = std_pitch_shift(3:end);
    std_pitch_shift_no_base(1) = 0;
    std_no_pitch_shift_no_base = std_no_pitch_shift(3:end);
    std_no_pitch_shift_no_base(1) = 0;
else
    mean_neg_pitch_shift_no_base = [0; mean_neg_pitch_shift];
    mean_pos_pitch_shift_no_base = [0; mean_pos_pitch_shift];
    std_neg_pitch_shift_no_base = [0 std_neg_pitch_shift];
    std_pos_pitch_shift_no_base = [0 std_pos_pitch_shift];
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 plots the change in pitch of lesioned shift birds both positive
%and negative shifts without baseline changes.

figure (1)
hold all
errorbar(0:14,mean_shift_pitch_shift_no_base,std_pitch_shift_no_base,'m','lineWidth',2,'DisplayName','-1 shift')
errorbar((0:14)+0.1,mean_no_pitch_shift_no_base,std_no_pitch_shift_no_base,'k','lineWidth',2,'DisplayName','+1 shift')
scatter(0:14,mean_shift_pitch_shift_no_base,30,'k','fill')
scatter((0:14)+0.1,mean_no_pitch_shift_no_base,30,'k','fill')
legend()


% 0 cents line
plot(xlim, [0 0], 'k','LineWidth',2);

ylabel('Pitch Shift from Baseline in semitones')
xlabel('Days of shift')

%%
%The following block of code is used to plot individual birds mean and
%error bars for lesioned and unlesioned birds.

clear
close all
clc

load('workspace_variables_unlesioned.mat') %Make sure to set the right path

num_birds = length(birds_to_use);
mean_neg_pitch_shift = zeros(num_base_days + num_shift_days,num_birds);
mean_pos_pitch_shift = zeros(num_base_days + num_shift_days,num_birds);
mean_no_pitch_shift = zeros(num_base_days + num_shift_days,num_birds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    for j = neg_shift_birds
        temp = [];
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
        else
            temp = vertcat(P(j).data{i,:});
        end
        mean_neg_pitch_shift(i,j) = mean(temp);
    end
    
    for j = pos_shift_birds
        temp = [];
        if lesioned_birds
            temp = horzcat(P(j).data{i,:});
        else
            temp = vertcat(P(j).data{i,:});
        end
        mean_pos_pitch_shift(i,j) = mean(temp);
    end
    
    if lesioned_birds
        for j = no_shift_birds
            temp = [];
            if lesioned_birds
                temp = horzcat(P(j).data{i,:});
            else
                temp = vertcat(P(j).data{i,:});
            end
            mean_no_pitch_shift(i,j) = mean(temp);
        end
    end
    
end

if lesioned_birds
    mean_neg_pitch_shift_no_base = mean_neg_pitch_shift(3:end,:);
    mean_neg_pitch_shift_no_base(1,:) = zeros(1,num_birds);
    mean_pos_pitch_shift_no_base = mean_pos_pitch_shift(3:end,:);
    mean_pos_pitch_shift_no_base(1,:) = zeros(1,num_birds);
    mean_no_pitch_shift_no_base = mean_no_pitch_shift(3:end,:);
    mean_no_pitch_shift_no_base(1,:) = zeros(1,num_birds);
%     std_neg_pitch_shift_no_base = std_neg_pitch_shift(3:end);
%     std_neg_pitch_shift_no_base(1) = 0;
%     std_pos_pitch_shift_no_base = std_pos_pitch_shift(3:end);
%     std_pos_pitch_shift_no_base(1) = 0;
%     std_no_pitch_shift_no_base = std_no_pitch_shift(3:end);
%     std_no_pitch_shift_no_base(1) = 0;
else
    mean_neg_pitch_shift_no_base = [zeros(1,num_birds); mean_neg_pitch_shift];
    mean_pos_pitch_shift_no_base = [zeros(1,num_birds); mean_pos_pitch_shift];
%     std_neg_pitch_shift_no_base = [0 std_neg_pitch_shift];
%     std_pos_pitch_shift_no_base = [0 std_pos_pitch_shift];
end

figure (1)
hold all
for j = neg_shift_birds
    plot(0:14,mean_neg_pitch_shift_no_base(:,j),'Color',[1 0.5 0.5],'LineWidth',2)
end
for j = pos_shift_birds
    plot(0:14,mean_pos_pitch_shift_no_base(:,j),'Color',[0.5 0.5 1],'LineWidth',2)
end
xlim([-0.2 14.2])
plot(xlim, [0 0], 'k','LineWidth',2);

ylabel('Pitch Shift from Baseline in semitones')
xlabel('Days of shift')

if lesioned_birds
    figure (2)
    hold all
    for j = no_shift_birds
        plot(0:14,mean_no_pitch_shift_no_base(:,j),'Color',[0.1*(j-7) 0.1*(j-7) 0.1*(j-7)],'LineWidth',2)
    end
    xlim([-0.2 14.2])
    plot(xlim, [0 0], 'k','LineWidth',2);

    ylabel('Pitch Shift from Baseline in semitones')
    xlabel('Days of shift')
end

%%
%The following block of code is for generating histograms of shift of
%individual syllables. So far, this plot has not been included in the paper
%but that might change based on revisions.

clear
close all
clc

load('workspace_variables_for_gordon.mat') %Make sure to set the right path

neg_shift_birds_syl_means = [];
pos_shift_birds_syl_means = [];
if lesioned_birds
    no_shift_birds_syl_means = [];
end

for i = neg_shift_birds
    num_syls = size(P(i).data,2);
    for j = 1:num_syls
        if lesioned_birds
            temp = horzcat(P(i).data{num_days-2:num_days,j});
            neg_shift_birds_syl_means = [neg_shift_birds_syl_means; mean(temp)];
        else
            temp = vertcat(P(i).data{num_days-2:num_days,j});
            neg_shift_birds_syl_means = [neg_shift_birds_syl_means; mean(temp)];
        end
    end
            
end

for i = pos_shift_birds
    num_syls = size(P(i).data,2);
    for j = 1:num_syls
        if lesioned_birds
            temp = horzcat(P(i).data{num_days-2:num_days,j});
            pos_shift_birds_syl_means = [pos_shift_birds_syl_means; mean(temp)];
        else
            temp = vertcat(P(i).data{num_days-2:num_days,j});
            pos_shift_birds_syl_means = [pos_shift_birds_syl_means; mean(temp)];
        end
    end
            
end

if lesioned_birds
    for i = no_shift_birds
        num_syls = size(P(i).data,2);
        for j = 1:num_syls
            if lesioned_birds
                temp = horzcat(P(i).data{num_days-2:num_days,j});
                no_shift_birds_syl_means = [no_shift_birds_syl_means; mean(temp)];
            end
        end

    end
end

bin_min = min([neg_shift_birds_syl_means; pos_shift_birds_syl_means]);
bin_max = max([neg_shift_birds_syl_means; pos_shift_birds_syl_means]);

bin_min = round(bin_min, 1) - 0.1;
bin_max = round(bin_max, 1) + 0.1;

bins_all = bin_min:0.1:bin_max;

neg_syls_hist = -1*histcounts(neg_shift_birds_syl_means,bins_all);
pos_syls_hist = histcounts(pos_shift_birds_syl_means,bins_all);

bins_mod = bins_all + 0.05;
bins_mod = bins_mod(1:end-1);

figure (1)
hold all
barh(bins_mod, neg_syls_hist,'r')
barh(bins_mod, pos_syls_hist,'b')
xlim([-5 5])
ylim([-2.5 1.5])
plot(xlim, [0 0], 'k--','LineWidth',2)

if lesioned_birds
    no_syls_hist = histcounts(no_shift_birds_syl_means,bins_all);

    figure (2)
    hold all
    barh(bins_mod, no_syls_hist,'k')
    plot(xlim, [0 0], 'k--','LineWidth',2)
end
