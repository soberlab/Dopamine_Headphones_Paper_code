%This is a script written to do stats on the 6-OHDA headphones results.

clear
close all
clc

%Will rename this file later: Has 6-OHDA shift and no shift birds data.
load('workspace_variables_for_gordon.mat')

rng(1); %For reproducibility of stats.

%Get summary statistics (N and means for days 12 through 14 of shift):
%The 4th element is for the shift birds combined (see next section).
%The 5th and 6th elements are for no shift group broken down by headphones
%vs no headphones (see 2 sections below).
means_summary = zeros(1,6);
means_n = zeros(1,6);
temp2_neg = [];
temp2_pos = [];
temp2_no = [];
for i = 15:num_days
    for j = neg_shift_birds
        temp = horzcat(P(j).data{i,:});
        temp2_neg = horzcat(temp2_neg,temp);
    end
    for j = pos_shift_birds
        temp = horzcat(P(j).data{i,:});
        temp2_pos = horzcat(temp2_pos,temp);
    end
    for j = no_shift_birds
        temp = horzcat(P(j).data{i,:});
        temp2_no = horzcat(temp2_no,temp);
    end
    
end
means_summary(1) = mean(temp2_neg);
means_summary(2) = mean(temp2_pos);
means_summary(3) = mean(temp2_no);
means_n(1) = numel(temp2_neg);
means_n(2) = numel(temp2_pos);
means_n(3) = numel(temp2_no);

%To calculate SEMs, we still need to do bootstrapping but now over all 3
%days.

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.

%num_days + 1 for the fnal column which takes the mean over days 12 - 14.
bootstrapping_matrix = zeros(nboot,num_days+1,max_syls,numel(birds_to_use));

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
            bootstrapping_matrix(n,num_days+1,j,i) = mean(bootstrapping_matrix(n,num_days-2:num_days,j,i));
            temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
            bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
        end
    end
end

%%
nboot2 = 10000; %Repeating 1000 times since this is the distribution used to calculate stats.
bootstats1 = zeros(nboot2,1);
bootstats2 = zeros(nboot2,1);
bootstats3 = zeros(nboot2,1);
for n = 1:nboot2
    
    temp_neg_birds = datasample(neg_shift_birds,length(neg_shift_birds));
    temp_data = [];
    for t = 1:length(temp_neg_birds)
        temp_syls = datasample(1:size(P(temp_neg_birds(t)).data,2),size(P(temp_neg_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_neg_birds(t)));
        end
    end
    bootstats1(n,1) = nanmean(temp_data);
    
    temp_pos_birds = datasample(pos_shift_birds,length(pos_shift_birds));
    temp_data = [];
    for t = 1:length(temp_pos_birds)
        temp_syls = datasample(1:size(P(temp_pos_birds(t)).data,2),size(P(temp_pos_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_pos_birds(t)));
        end
    end
    bootstats2(n,1) = nanmean(temp_data);
    
    temp_no_birds = datasample(no_shift_birds,length(no_shift_birds));
    temp_data = [];
    for t = 1:length(temp_no_birds)
        temp_syls = datasample(1:size(P(temp_no_birds(t)).data,2),size(P(temp_no_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_no_birds(t)));
        end
    end
    bootstats3(n,1) = nanmean(temp_data);
end

means_sem = zeros(1,6);
means_sem(1) = nanstd(bootstats1);
means_sem(2) = nanstd(bootstats2);
means_sem(3) = nanstd(bootstats3);

%%
%In the below block of code, we repeat the analysis as above but now we are
%interested in simply shift versus no shift birds.

temp2_shift = [];
shift_birds = [neg_shift_birds pos_shift_birds];
for i = 15:num_days
    for j = shift_birds
        temp = horzcat(P(j).data{i,:});
        temp2_shift = horzcat(temp2_shift,temp);
    end
    
end
means_summary(4) = mean(temp2_shift);
means_n(4) = numel(temp2_shift);

%To calculate SEMs, we still need to do bootstrapping but now over all 3
%days. We are going to use the same bootstrapping_matrix from section 1.

nboot2 = 10000; 
bootstats4 = zeros(nboot2,1);
for n = 1:nboot2
    
    temp_shift_birds = datasample(shift_birds,length(shift_birds));
    temp_data = [];
    for t = 1:length(temp_shift_birds)
        temp_syls = datasample(1:size(P(temp_shift_birds(t)).data,2),size(P(temp_shift_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_shift_birds(t)));
        end
    end
    bootstats4(n,1) = nanmean(temp_data);
    
end

means_sem(4) = nanstd(bootstats4);

%%
%In the below block of code, we repeat the analysis as above but now we are
%interested in no shift birds broken down by headphones vs not.

no_shift_no_headp = 8:10;
no_shift_headp = 11:15;

temp2_no_headp = [];
temp2_no_no_headp = [];
for i = 15:num_days
    for j = no_shift_headp
        temp = horzcat(P(j).data{i,:});
        temp2_no_headp = horzcat(temp2_no_headp,temp);
    end
    for j = no_shift_no_headp
        temp = horzcat(P(j).data{i,:});
        temp2_no_no_headp = horzcat(temp2_no_no_headp,temp);
    end
    
end
means_summary(5) = mean(temp2_no_headp);
means_summary(6) = mean(temp2_no_no_headp);
means_n(5) = numel(temp2_no_headp);
means_n(6) = numel(temp2_no_no_headp);

nboot2 = 10000; 
bootstats5 = zeros(nboot2,1);
bootstats6 = zeros(nboot2,1);
for n = 1:nboot2
    
    temp_shift_birds = datasample(no_shift_headp,length(no_shift_headp));
    temp_data = [];
    for t = 1:length(temp_shift_birds)
        temp_syls = datasample(1:size(P(temp_shift_birds(t)).data,2),size(P(temp_shift_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_shift_birds(t)));
        end
    end
    bootstats5(n,1) = nanmean(temp_data);
    
    temp_shift_birds = datasample(no_shift_no_headp,length(no_shift_no_headp));
    temp_data = [];
    for t = 1:length(temp_shift_birds)
        temp_syls = datasample(1:size(P(temp_shift_birds(t)).data,2),size(P(temp_shift_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_shift_birds(t)));
        end
    end
    bootstats6(n,1) = nanmean(temp_data);
    
end

means_sem(5) = nanstd(bootstats5);
means_sem(6) = nanstd(bootstats6);

%%
%The following section calculates the statistical probabilities on the
%above groups based on the bootstrapped samples.
%First, we calculate the probability of the no shift birds being less than
%zero.

p_no_shift_zero = sum(bootstats3 >= 0)/length(bootstats3);

%The more involved part is to calculate the joint probability distribution
%between two groups and see if they are different. We want to do this
%between the populations of no shift birds - headphones vs no headphones.
%But this serves as a template for other comparisons also. 

%Done and converted into a function - get_direct_prob.m.

%So, calculate between no shift headphones and no headphones:
p_headp_vs_no_headp = get_direct_prob(bootstats5,bootstats6);
%The result, the prob of no headphones being greater than or equal to
%headphones = 0.91

%Calculate between neg shift and pos shift:
[p_neg_shift_vs_pos_shift, p_dist] = get_direct_prob(bootstats1, bootstats2);
%The result, the prob of pos_shift being greater than or equal to neg_shift
%= 0.26 (N = 8)
%The result, the prob of pos_shift being greater than or equal to neg_shift
%= 0.44 (N = 7)

%Calculate between shift and no shift:
p_shift_vs_no_shift = get_direct_prob(bootstats4, bootstats3);
%The result, the prob of no_shift being greater than or equal to shift =
%0.79 (N = 8)
%The result, the prob of no_shift being greater than or equal to shift =
%0.91 (N = 7)

%Calculate neg shift vs no shift and pos shift vs no shift individually:
p_neg_vs_no = get_direct_prob(bootstats1, bootstats3);
%The result, the prob of no_shift being greater than or equal to neg shift =
%0.62 (N = 8)
%The result, the prob of no_shift being greater than or equal to neg shift =
%0.80 (N = 7)

p_pos_vs_no = get_direct_prob(bootstats2, bootstats3);
%The result, the prob of no_shift being greater than or equal to pos shift =
%0.91 (N = 8) Same for (N = 7)



%%
%The following section is for unlesioned birds and analysis.
clear
close all
clc

load('workspace_variables_unlesioned.mat')

rng(1); %For reproducibility

means_summary = zeros(1,2);
means_n = zeros(1,2);
temp2_neg = [];
temp2_pos = [];
for i = 12:num_days
    for j = neg_shift_birds
        temp = vertcat(P(j).data{i,:});
        temp2_neg = vertcat(temp2_neg,temp);
    end
    for j = pos_shift_birds
        temp = vertcat(P(j).data{i,:});
        temp2_pos = vertcat(temp2_pos,temp);
    end
    
end
means_summary(1) = mean(temp2_neg);
means_summary(2) = mean(temp2_pos);
means_n(1) = numel(temp2_neg);
means_n(2) = numel(temp2_pos);

%To calculate SEMs, we still need to do bootstrapping but now over all 3
%days.

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.

%num_days + 1 for the final column which takes the mean over days 12 - 14.
bootstrapping_matrix = zeros(nboot,num_days+1,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(P(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days      %Over all days
            temp = P(i).data{k,j};
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
        for n = 1:nboot
            bootstrapping_matrix(n,num_days+1,j,i) = mean(bootstrapping_matrix(n,num_days-2:num_days,j,i));
        end
    end
end

nboot2 = 10000;
bootstats1 = zeros(nboot2,1);
bootstats2 = zeros(nboot2,1);
for n = 1:nboot2
    
    temp_neg_birds = datasample(neg_shift_birds,length(neg_shift_birds));
    temp_data = [];
    for t = 1:length(temp_neg_birds)
        temp_syls = datasample(1:size(P(temp_neg_birds(t)).data,2),size(P(temp_neg_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_neg_birds(t)));
        end
    end
    bootstats1(n,1) = nanmean(temp_data);
    
    temp_pos_birds = datasample(pos_shift_birds,length(pos_shift_birds));
    temp_data = [];
    for t = 1:length(temp_pos_birds)
        temp_syls = datasample(1:size(P(temp_pos_birds(t)).data,2),size(P(temp_pos_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_pos_birds(t)));
        end
    end
    bootstats2(n,1) = nanmean(temp_data);
    
end

means_sem = zeros(1,2);
means_sem(1) = nanstd(bootstats1);
means_sem(2) = nanstd(bootstats2);

p_pos_shift_vs_neg_shift = get_direct_prob(bootstats1, bootstats2);

%%
%The following section is for unlesioned birds washout analysis.
clear
close all
clc

load('workspace_variables_unlesioned_washout.mat')
rng(1);

means_summary = zeros(1,2);
means_n = zeros(1,2);
temp2_neg = [];
temp2_pos = [];
for i = 6:num_days
    for j = neg_shift_birds
        temp = vertcat(Q(j).data{i,:});
        temp2_neg = vertcat(temp2_neg,temp);
    end
    for j = pos_shift_birds
        temp = vertcat(Q(j).data{i,:});
        temp2_pos = vertcat(temp2_pos,temp);
    end
    
end
means_summary(1) = mean(temp2_neg);
means_summary(2) = mean(temp2_pos);
means_n(1) = numel(temp2_neg);
means_n(2) = numel(temp2_pos);

%To calculate SEMs, we still need to do bootstrapping but now over all 3
%days.

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.

%num_days + 1 for the final column which takes the mean over days 12 - 14.
bootstrapping_matrix = zeros(nboot,num_days+1,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(Q(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days      %Over all days
            temp = Q(i).data{k,j};
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            for n = 1:nboot     %For number of times to bootstrap
                bootstrapping_matrix(n,k,j,i) = datasample(temp,1);
            end
        end
        for n = 1:nboot
            bootstrapping_matrix(n,num_days+1,j,i) = mean(bootstrapping_matrix(n,num_days-1:num_days,j,i));
        end
    end
end

nboot2 = 10000; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,1);
bootstats2 = zeros(nboot2,1);
for n = 1:nboot2
    
    temp_neg_birds = datasample(neg_shift_birds,length(neg_shift_birds));
    temp_data = [];
    for t = 1:length(temp_neg_birds)
        temp_syls = datasample(1:size(Q(temp_neg_birds(t)).data,2),size(Q(temp_neg_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_neg_birds(t)));
        end
    end
    bootstats1(n,1) = nanmean(temp_data);
    
    temp_pos_birds = datasample(pos_shift_birds,length(pos_shift_birds));
    temp_data = [];
    for t = 1:length(temp_pos_birds)
        temp_syls = datasample(1:size(Q(temp_pos_birds(t)).data,2),size(Q(temp_pos_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_pos_birds(t)));
        end
    end
    bootstats2(n,1) = nanmean(temp_data);
    
end

means_sem = zeros(1,2);
means_sem(1) = nanstd(bootstats1);
means_sem(2) = nanstd(bootstats2);

p_pos_shift_zero = sum(bootstats2 > 0)/length(bootstats2);
p_neg_shift_zero = sum(bootstats1 < 0)/length(bootstats1);


%%
%The following section is for lesioned birds washout analysis.

clear
close all
clc

%Will rename this file later: Has 6-OHDA shift and no shift birds data.
load('workspace_variables_lesioned_washout.mat')
rng(1);

%Get summary statistics (N and means for days 6 and 7 of washout):
means_summary = zeros(1,2);
means_n = zeros(1,2);
temp2_neg = [];
temp2_pos = [];
for i = (num_days-1):num_days
    for j = neg_shift_birds
        temp = horzcat(Q(j).data{i,:});
        temp2_neg = horzcat(temp2_neg,temp);
    end
    for j = pos_shift_birds
        temp = horzcat(Q(j).data{i,:});
        temp2_pos = horzcat(temp2_pos,temp);
    end
    
end
means_summary(1) = mean(temp2_neg);
means_summary(2) = mean(temp2_pos);
means_n(1) = numel(temp2_neg);
means_n(2) = numel(temp2_pos);

%To calculate SEMs, we still need to do bootstrapping but now over all 3
%days.

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.

%num_days + 1 for the fnal column which takes the mean over days 12 - 14.
bootstrapping_matrix = zeros(nboot,25,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(Q(i).data,2) %Over number of syllables for each bird
        for k = 1:size(Q(i).data_hz,1)      %Over all days
            temp = Q(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
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
            bootstrapping_matrix(n,end_of_shift:25,j,i) = bootstrapping_matrix(n,end_of_shift:25,j,i) - bootstrapping_matrix(n,end_of_shift,j,i);
            bootstrapping_matrix(n,25,j,i) = mean(bootstrapping_matrix(n,23:24,j,i));
        end
    end
end

nboot2 = 10000; %Only repeating 300 times since the actual variance is in the previous step.
bootstats1 = zeros(nboot2,1);
bootstats2 = zeros(nboot2,1);
for n = 1:nboot2
    
    temp_neg_birds = datasample(neg_shift_birds,length(neg_shift_birds));
    temp_data = [];
    for t = 1:length(temp_neg_birds)
        temp_syls = datasample(1:size(Q(temp_neg_birds(t)).data,2),size(Q(temp_neg_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_neg_birds(t)));
        end
    end
    bootstats1(n,1) = nanmean(temp_data);
    
    temp_pos_birds = datasample(pos_shift_birds,length(pos_shift_birds));
    temp_data = [];
    for t = 1:length(temp_pos_birds)
        temp_syls = datasample(1:size(Q(temp_pos_birds(t)).data,2),size(Q(temp_pos_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_pos_birds(t)));
        end
    end
    bootstats2(n,1) = nanmean(temp_data);
    
end

means_sem = zeros(1,2);
means_sem(1) = nanstd(bootstats1);
means_sem(2) = nanstd(bootstats2);

p_pos_shift_zero = sum(bootstats2 > 0)/length(bootstats2);
p_neg_shift_zero = sum(bootstats1 < 0)/length(bootstats1);

%%
%We need a section to perform a check on lesioned headphones no shift birds
%for shift period versus washout period.

clc
clear
close all

lesioned_birds = 1;

if lesioned_birds
    birds_to_use = 35:39; %Only 6-OHDA headphones no shift birds
    no_shift_birds = 1:5;
    %Set these:
    num_base_days = 3;
    num_shift_days = 14;
    num_washout_days = 2;
    num_days = num_base_days + num_shift_days + num_washout_days;
    max_syls = 12;
end

counter = 1;
for i = birds_to_use
    P(counter) = load_bird_params_2018(i,'washout_no_shift');
    counter = counter + 1;
end

rng(1); %For reproducibility of stats.

%We only need data from Days 15-17 and Days 18-19 (correspond to days 6 and
%7 in washout respectively).
means_summary = zeros(1,2);
temp2_no = [];
temp2_no_washout = [];
for i = 15:17
    for j = no_shift_birds
        temp = horzcat(P(j).data{i,:});
        temp2_no = horzcat(temp2_no,temp);
    end
    
end
means_summary(1) = mean(temp2_no);

for i = 18:19
    for j = no_shift_birds
        temp = horzcat(P(j).data{i,:});
        temp2_no_washout = horzcat(temp2_no_washout,temp);
    end
    
end
means_summary(2) = mean(temp2_no_washout);

%To calculate SEMs, we still need to do bootstrapping but now over all 3
%days.

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.

%num_days + 2 for the final column which takes the mean over days 12 - 14
%and for last 2 days of washout.
bootstrapping_matrix = zeros(nboot,num_days+2,max_syls,numel(birds_to_use));

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
            bootstrapping_matrix(n,num_days+1,j,i) = mean(bootstrapping_matrix(n,num_days-4:num_days-2,j,i));
            bootstrapping_matrix(n,num_days+2,j,i) = mean(bootstrapping_matrix(n,num_days-1:num_days,j,i));
            temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
            bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
        end
    end
end

nboot2 = 10000; %Repeating 1000 times since this is the distribution used to calculate stats.
bootstats1 = zeros(nboot2,1);
bootstats2 = zeros(nboot2,1);
for n = 1:nboot2
    
    temp_no_birds = datasample(no_shift_birds,length(no_shift_birds));
    temp_data = [];
    for t = 1:length(temp_no_birds)
        temp_syls = datasample(1:size(P(temp_no_birds(t)).data,2),size(P(temp_no_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end-1,temp_syls(s),temp_no_birds(t)));
        end
    end
    bootstats1(n,1) = nanmean(temp_data);
    
    temp_no_birds = datasample(no_shift_birds,length(no_shift_birds));
    temp_data = [];
    for t = 1:length(temp_no_birds)
        temp_syls = datasample(1:size(P(temp_no_birds(t)).data,2),size(P(temp_no_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_no_birds(t)));
        end
    end
    bootstats2(n,1) = nanmean(temp_data);
    
end

%Calculate neg shift vs no shift and pos shift vs no shift individually:
p_shift_vs_washout = get_direct_prob(bootstats2, bootstats1);
%The result, the prob of washout being greater than or equal to shift =
%0.66

%%
%This section is to check for no shift birds again with subtracted washout
%values compared to zero.
clc
clear
close all

lesioned_birds = 1;

if lesioned_birds
    birds_to_use = 35:39; %Only 6-OHDA headphones no shift birds
    no_shift_birds = 1:5;
    %Set these:
    num_base_days = 3;
    num_shift_days = 14;
    num_washout_days = 7;
    num_days = 8;
    num_bootstrapped_days = num_base_days + num_shift_days + num_washout_days;
    max_syls = 12;
    end_of_shift = 17;
end

counter = 1;
for i = birds_to_use
    P(counter) = load_bird_params_2018(i,'washout');
    counter = counter + 1;
end

rng(1); %For reproducibility of stats.

%We only need data from Days 6 and 7 which correspond to columns 7 and 8
%respectively.
means_summary = zeros(1);
temp2_no = [];
temp2_no_washout = [];
for i = 7:8
    for j = no_shift_birds
        temp = horzcat(P(j).data{i,:});
        temp2_no = horzcat(temp2_no,temp);
    end
    
end
means_summary(1) = mean(temp2_no);

%To calculate SEMs, we still need to do bootstrapping but now over all 3
%days.

%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.

%num_days + 2 for the final column which takes the mean over days 12 - 14
%and for last 2 days of washout.
bootstrapping_matrix = zeros(nboot,num_bootstrapped_days+1,max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(P(i).data,2) %Over number of syllables for each bird
        for k = 1:num_bootstrapped_days      %Over all days
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
            bootstrapping_matrix(n,end_of_shift:25,j,i) = bootstrapping_matrix(n,end_of_shift:25,j,i) - bootstrapping_matrix(n,end_of_shift,j,i);
            bootstrapping_matrix(n,25,j,i) = mean(bootstrapping_matrix(n,23:24,j,i));
        end
    end
end

nboot2 = 10000; %Repeating 1000 times since this is the distribution used to calculate stats.
bootstats1 = zeros(nboot2,1);
for n = 1:nboot2
    
    temp_no_birds = datasample(no_shift_birds,length(no_shift_birds));
    temp_data = [];
    for t = 1:length(temp_no_birds)
        temp_syls = datasample(1:size(P(temp_no_birds(t)).data,2),size(P(temp_no_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_no_birds(t)));
        end
    end
    bootstats1(n,1) = nanmean(temp_data);
    
end

%Calculate neg shift vs no shift and pos shift vs no shift individually:
sem_summary = std(bootstats1);
p_no_shift_zero = sum(bootstats1 >= 0)/length(bootstats1);
%The probability of this being greater than or equal to zero is 0.2168.

%%
%We need a section to compare last days of shift versus last days of
%washout directly for unlesioned and lesioned shift groups.

clear
close all
clc

%However all reported stats previously have been performed on subtracted
%data for washout. Therefore, we have to reload the data and compute new
%bootstrapped samples (we can keep the bootstrapped samples from the shift
%period though).

%So first, we are going to compute and same the equivalent bootstats
%matrices for lesioned and unlesioned birds since it requires the old
%fashioned loading of raw data.

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

if lesioned_birds
    mean_neg_pitch_shift = zeros(2);
    mean_pos_pitch_shift = zeros(2);

    temp2_neg = [];
    temp2_pos = [];
    for i = (end_of_shift-2):end_of_shift

        for j = neg_shift_birds
            if lesioned_birds
                temp = horzcat(Q(j).data{i,:});
                temp2_neg = horzcat(temp2_neg,temp);
            else
                temp = vertcat(Q(j).data{i,:});
                temp2_neg = vertcat(temp2_neg,temp);
            end
        end

        for j = pos_shift_birds
            if lesioned_birds
                temp = horzcat(Q(j).data{i,:});
                temp2_pos = horzcat(temp2_pos,temp);
            else
                temp = vertcat(Q(j).data{i,:});
                temp2_pos = vertcat(temp2_pos,temp);
            end
        end

    end
    mean_neg_pitch_shift(1) = mean(temp2_neg);
    mean_pos_pitch_shift(1) = mean(temp2_pos);
end

temp2_neg = [];
temp2_pos = [];
for i = (num_days-1):num_days
    
    for j = neg_shift_birds
        if lesioned_birds
            temp = horzcat(Q(j).data{i,:});
            temp2_neg = horzcat(temp2_neg,temp);
        else
            temp = vertcat(Q(j).data{i,:});
            temp2_neg = vertcat(temp2_neg,temp);
        end
    end
    
    for j = pos_shift_birds
        if lesioned_birds
            temp = horzcat(Q(j).data{i,:});
            temp2_pos = horzcat(temp2_pos,temp);
        else
            temp = vertcat(Q(j).data{i,:});
            temp2_pos = vertcat(temp2_pos,temp);
        end
    end
   
end
if lesioned_birds
    mean_neg_pitch_shift(2) = mean(temp2_neg);
    mean_pos_pitch_shift(2) = mean(temp2_pos);
else
    mean_neg_pitch_shift = mean(temp2_neg);
    mean_pos_pitch_shift = mean(temp2_pos);
end
%Bootstrapping part
nboot = 1000; %No of times to resample for bootstrapping.
if lesioned_birds
    num_days_bootstrap = size(Q(1).data_hz,1);
    max_syls = 12;
else
    num_days_bootstrap = size(Q(1).data,1);
    max_syls = 9;
end

if lesioned_birds
    bootstrapping_matrix = zeros(nboot,(num_days_bootstrap+2),max_syls,numel(birds_to_use));
else
    bootstrapping_matrix = zeros(nboot,(num_days_bootstrap+1),max_syls,numel(birds_to_use));
end

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
            bootstrapping_matrix(n,num_days_bootstrap+1,j,i) = mean(bootstrapping_matrix(n,num_days_bootstrap-1:num_days_bootstrap,j,i));
            if lesioned_birds
                bootstrapping_matrix(n,num_days_bootstrap+2,j,i) = mean(bootstrapping_matrix(n,end_of_shift-2:end_of_shift,j,i));
                temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
                bootstrapping_matrix(n,:,j,i) = 12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            end
            
        end
    end
end


bootstrapping_matrix = bootstrapping_matrix(:,end_of_shift:end,:,:);

%Now do the actual bootstrapping: 
nboot2 = 10000; 
bootstats1_s = zeros(nboot2,1);
bootstats2_s = zeros(nboot2,1);
bootstats1_w = zeros(nboot2,1);
bootstats2_w = zeros(nboot2,1);
if lesioned_birds
    std_neg_pitch_shift = zeros(2);
    std_pos_pitch_shift = zeros(2);
end
for n = 1:nboot2
    
    temp_neg_birds = datasample(neg_shift_birds,length(neg_shift_birds));
    temp_data = [];
    temp_data2 = [];
    for t = 1:length(temp_neg_birds)
        temp_syls = datasample(1:size(Q(temp_neg_birds(t)).data,2),size(Q(temp_neg_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_neg_birds(t)));
            temp_data2 = vertcat(temp_data2,bootstrapping_matrix(temp_pulls,(end-1),temp_syls(s),temp_neg_birds(t)));
        end
    end
    bootstats1_s(n,1) = nanmean(temp_data,1);
    bootstats1_w(n,1) = nanmean(temp_data2,1);
    
    temp_pos_birds = datasample(pos_shift_birds,length(pos_shift_birds));
    temp_data = [];
    temp_data2 = [];
    for t = 1:length(temp_pos_birds)
        temp_syls = datasample(1:size(Q(temp_pos_birds(t)).data,2),size(Q(temp_pos_birds(t)).data,2));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_pos_birds(t)));
            temp_data2 = vertcat(temp_data2,bootstrapping_matrix(temp_pulls,(end-1),temp_syls(s),temp_pos_birds(t)));
        end
    end
    bootstats2_s(n,1) = nanmean(temp_data,1);
    bootstats2_w(n,1) = nanmean(temp_data2,1);
end
if lesioned_birds
    std_neg_pitch_shift(1) = nanstd(bootstats1_s);
    std_pos_pitch_shift(1) = nanstd(bootstats2_s);
    std_neg_pitch_shift(2) = nanstd(bootstats1_w);
    std_pos_pitch_shift(2) = nanstd(bootstats2_w);
else
    std_neg_pitch_shift = nanstd(bootstats1_s);
    std_pos_pitch_shift = nanstd(bootstats2_s);
end

% p1 = get_direct_prob(lesioned_shift.bootstats1, lesioned_washout.bootstats1)
% p2 = get_direct_prob(lesioned_shift.bootstats2, lesioned_washout.bootstats2)
% p3 = get_direct_prob(unlesioned_shift.bootstats1, unlesioned_washout.bootstats1)
% p4 = get_direct_prob(unlesioned_shift.bootstats2, unlesioned_washout.bootstats2)

%%
%The following section is used for Linear Mixed Model analysis of lesioned
%shift data.

clear
close all
clc

%Will rename this file later: Has 6-OHDA shift and no shift birds data.
load('workspace_variables_for_gordon_n_7.mat')

%Get summary statistics (N and means for days 12 through 14 of shift):
%The 4th element is for the shift birds combined (see next section).
%The 5th and 6th elements are for no shift group broken down by headphones
%vs no headphones (see 2 sections below).

temp_table = [];

for i = 15:num_days
    for j = neg_shift_birds
        for k = 1:size(P(j).data,2)
            temp = P(j).data{i,k}';
            temp_bird = j*ones(length(temp),1);
            temp_syl = k*ones(length(temp),1);
            temp_cond = -1*ones(length(temp),1);
            temp_cond_2 = 1*ones(length(temp),1);
            temp_table = vertcat(temp_table, [temp temp_bird temp_syl temp_cond temp_cond_2]);
        end
    end
    for j = pos_shift_birds
        for k = 1:size(P(j).data,2)
            temp = P(j).data{i,k}';
            temp_bird = j*ones(length(temp),1);
            temp_syl = k*ones(length(temp),1);
            temp_cond = 1*ones(length(temp),1);
            temp_table = vertcat(temp_table, [temp temp_bird temp_syl temp_cond temp_cond]);
        end
    end
    for j = no_shift_birds
        for k = 1:size(P(j).data,2)
            temp = P(j).data{i,k}';
            temp_bird = j*ones(length(temp),1);
            temp_syl = k*ones(length(temp),1);
            temp_cond = -1*ones(length(temp),1);
            temp_cond_2 = zeros(length(temp),1);
            temp_table = vertcat(temp_table, [temp temp_bird temp_syl temp_cond temp_cond_2]);
        end
    end
    
end

tbl = table(temp_table(:,1),temp_table(:,2), temp_table(:,3), temp_table(:,4), temp_table(:,5),...
    'VariableNames',{'Pitch','Bird','Syllable','Shift_cond','Shift'});

%The linear model we fit will have the shift condition as the fixed effect
%and the bird and syllable as variable effects that are hierarchical.

lme1 = fitlme(tbl,'Pitch ~ Shift_cond + (1|Bird) + (1|Bird:Syllable)');

%Also fit the model to distinguish between shift versus no shift using the
%Shift variable defined above.

lme2 = fitlme(tbl,'Pitch ~ Shift + (1|Bird) + (1|Bird:Syllable)');

%%
%The following section is used for Linear Mixed Model analysis of unlesioned
%shift data.

clear
close all
clc

%Will rename this file later: Has 6-OHDA shift and no shift birds data.
load('workspace_variables_unlesioned.mat')

%Get summary statistics (N and means for days 12 through 14 of shift):
%The 4th element is for the shift birds combined (see next section).
%The 5th and 6th elements are for no shift group broken down by headphones
%vs no headphones (see 2 sections below).

temp_table = [];

for i = 12:num_days
    for j = neg_shift_birds
        for k = 1:size(P(j).data,2)
            temp = P(j).data{i,k};
            temp_bird = j*ones(length(temp),1);
            temp_syl = k*ones(length(temp),1);
            temp_cond = -1*ones(length(temp),1);
            if size(temp,1)==1
                temp_table = vertcat(temp_table, [temp' temp_bird temp_syl temp_cond]);
            else
                temp_table = vertcat(temp_table, [temp temp_bird temp_syl temp_cond]);
            end
        end
    end
    for j = pos_shift_birds
        for k = 1:size(P(j).data,2)
            temp = P(j).data{i,k}';
            temp_bird = j*ones(length(temp),1);
            temp_syl = k*ones(length(temp),1);
            temp_cond = 1*ones(length(temp),1);
            if size(temp,1)==1
                temp_table = vertcat(temp_table, [temp' temp_bird temp_syl temp_cond]);
            else
                temp_table = vertcat(temp_table, [temp temp_bird temp_syl temp_cond]);
            end
        end
    end
    
end

tbl = table(temp_table(:,1),temp_table(:,2), temp_table(:,3), temp_table(:,4),...
    'VariableNames',{'Pitch','Bird','Syllable','Shift_cond'});

%The linear model we fit will have the shift condition as the fixed effect
%and the bird and syllable as variable effects that are hierarchical.

lme1 = fitlme(tbl,'Pitch ~ Shift_cond + (1|Bird) + (1|Bird:Syllable)');


