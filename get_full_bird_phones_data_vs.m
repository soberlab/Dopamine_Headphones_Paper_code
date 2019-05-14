function [ pitches ] = get_full_bird_phones_data_vs(birdname, num_days, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin > 2
    cond = varargin{1};
end
switch birdname
    case 'bl20gr152'
        syls = {'w';'e';'a';'s';'d'};
        n_syl_in_sequence = [1 1 1 1 1];
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_no_headphones_birds\bl20gr152_labeled_song';
        shift_dir = '';
        washout_dir = '';
        
    case 'pu63pu64'
        syls = {'w';'e';'a';'s';'d'};
        n_syl_in_sequence = [1 1 1 1 1];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_no_headphones_birds\pu63pu64_labeled_song\pu63pu64_labeled_baseline';
        shift_dir = '';
        washout_dir = '';
        
    case 'bl23lb23'
        syls = {'a';'s';'d';'f';'r';'w'};
        n_syl_in_sequence = [1 1 1 1 1 1];
        
        baseline_dir = 'C:\DATA\bl23lb23_data\bl23lb23_labeled_song\bl23lb23_post6ohda_baseline';
        shift_dir = '';
        washout_dir = '';
        
    case 'bl36yw96'
        syls = {'s';'d';'e';'f';'j';'l'};
        n_syl_in_sequence = [1 1 1 1 1 1];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\bl36yw96\bl36yw96_baseline_labeled_songs';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\bl36yw96\bl36yw96_neg100shift_labeled_songs';
        washout_dir = '';
        
    case 'bk9gy9'
        syls = {'a';'s';'d';'j';'k';'l'};
        n_syl_in_sequence = [1 1 1 1 1 1];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\bk9gy9\bk9gy9_labeled_baseline';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\bk9gy9\bk9gy9_labeled_neg100shift';
        washout_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\bk9gy9\bk9gy9_labeled_washout';
%         washout_dir = '';
        
    case 'bl47yw80'
        syls = {'a';'s';'d';'f';'j';'k';'l';'i'};
        n_syl_in_sequence = [1 1 1 1 1 1 1 1];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\bl47yw80\bl47yw80_baseline_labeled_songs'; % original files stored on other computer but these backup files are identical. you must have backup drive plugged in
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\bl47yw80\bl47yw80_100shift_labeled_songs';    % original files stored on other computer but these backup files are identical. you must have backup drive plugged in
        washout_dir = '';
        
    case 'lb160rd190'
        syls = {'s';'d';'f';'e';'h';'j';'k';'t';'l';'i';'a';'g'};
        %syls = {'f';'e';'l';'a';'g'}; % only plot subset of syls.
        n_syl_in_sequence = [1 1 1 1 1 1 1 1 1 1 1 1];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\lb160rd190\lb160rd190_baseline_labeled_songs'; % original files stored on other computer but these backup files are identical. you must have backup drive plugged in
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\lb160rd190\lb160rd190_neg100shift_labeled_songs'; % original files stored on other computer but these backup files are identical. you must have backup drive plugged in
        washout_dir = '';
       
    case 'lb199gy78'
        syls = {'a';'s';'d';'f';'i'};
        n_syl_in_sequence = [1 1 1 1 1];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\lb199gy78\lb199gy78_baseline_labeled_songs';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\lb199gy78\lb199gy78_100shift_labeled_songs';
        washout_dir = '';
        
    case 'pu81pu82'
        syls = {'a';'s';'d';'f';'g'};
        n_syl_in_sequence = [1 1 1 1 1 ];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\pu81pu82\pu81pu82_labeled_song_baseline';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\pu81pu82\pu81pu82_labeled_song_baseline_drift';
%         washout_dir = '';
        washout_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\pu81pu82\pu81pu82_labeled_song_washout';
        
    case 'bk21bk22'
        all_syls = 0;
        if all_syls
            syls = {'a';'s';'d';'e';'r';'f';'g';'w'};
            n_syl_in_sequence = [1 1 1 1 1 1 1 1];
            colors = 'rrgbbmck';
            linewidths = [2 3 2 2 3 2 2 2];
        else
            syls = {'a';'s';'e';'r'};
            n_syl_in_sequence = [1 1 1 1];
            
        end
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\bk21bk22\bk21bk22_labeled_songs_baseline';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\bk21bk22\bk21bk22_labeled_songs_baseline_drift';
%         washout_dir = '';
        washout_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\bk21bk22\bk21bk22_labeled_washout';
       
    case 'br41br42'
        syls = {'a';'s';'d';'q';'f'};
        n_syl_in_sequence = [1 1 1 1 1];
        
        baseline_dir = 'C:\DATA\br41br42_data\br41br42_SHAM_neg100shift_labeled_songs\br41br42_baseline_labeled_songs';
        shift_dir = 'C:\DATA\br41br42_data\br41br42_SHAM_neg100shift_labeled_songs\br41br42_neg100shift_labeled_songs';
        washout_dir = '';
        
    case 'lb140yw4'
        syls = {'a';'s';'d';'f';'g'};
        n_syl_in_sequence = [1 1 1 1 1];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\lb140yw4\lb140yw4_baseline_labeled_songs';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\lb140yw4\lb140yw4_baseline_drift_labeled_songs';
%         washout_dir = '';
        washout_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\lb140yw4\lb140yw4_washout';
        
    case 'bk24gy24'
        syls = {'a';'s';'d';'k';'q';'w'};
        n_syl_in_sequence = [1 1 1 1 1 1];
       
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\bk24gy24\bk24gy24_labeled_songs_baseline';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\bk24gy24\bk24gy24_labeled_baseline_drift';
%         washout_dir = '';
        washout_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\bk24gy24\bk24gy24_labeled_washout';
        
    case 'lb24or124'
        syls = {'a';'s';'d';'f';'q'};
        n_syl_in_sequence = [1 1 1 1 1];
       
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\lb24or124\lb24or124_labeled_songs_baseline';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\lb24or124\lb24or124_labeled_songs_baseline_drift';
%         washout_dir = '';
        washout_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_noshift_birds\lb24or124\lb24or124_labeled_washout';
        
    case 'lb152rd152'
        if 0 % if plotting all syls that I labeled
            syls = {'a';'s';'d';'j';'k';'l';'f';'i'};
            n_syl_in_sequence = [1 1 1 1 1 1 1 1];
            colors = ['rrrgbbmc']; % similar-looking syls (asd, kl) plotted in the same color
            linewidths = [2 3 4 2 2 3 2 2];
            baseline_dir = 'H:\lb152rd152_6ohda_lesion_project\lb152rd152_labeled_songs\lb152rd152_baseline_labeled_songs';
            shift_dir = 'H:\lb152rd152_6ohda_lesion_project\lb152rd152_labeled_songs\lb152rd152_neg100shift_labeled_songs';
        else %if plotting q = asd combined and w = kl combined instead of separate a, s, d, k, l.
            %q = a, s, d combined into one syllable. It can be hard to tell difference between the three so I am making a 2nd set of plots of these as if they were one syllable.
            %w = k, l combined into one syllable. It can be hard to tell difference between the two so I am also making a 2nd set of plots of these as if they were one syllable.
            syls = {'q';'w';'j';'f';'i'};
            n_syl_in_sequence = [1 1 1 1 1];
           
            baseline_dir = 'H:\lb152rd152_6ohda_lesion_project\lb152rd152_labeled_songs_COMBINED_SYLS\lb152rd152_baseline_labeled_songs'; % Note the different folders used! One has asd/kl, one has q/w
            shift_dir = 'H:\lb152rd152_6ohda_lesion_project\lb152rd152_labeled_songs_COMBINED_SYLS\lb152rd152_neg100shift_labeled_songs'; % Note the different folders used! One has asd/kl, one has q/w
        end
        washout_dir = '';
        
    case 'gy46pu6'
        syls = {'a';'s';'f';'g'};
        n_syl_in_sequence = [1 1 1 1];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\gy46pu6\gy46pu6_labeled_baseline';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\gy46pu6\gy46pu6_labeled_pos100_shift';
        washout_dir = 'F:\DATA\gy46pu6_data\gy46pu6_labeled_songs\gy46pu6_labeled_washout';
%         washout_dir = '';
        
    case 'bl26bl27'
        syls = {'a';'s';'d';'f';'q';'w'};
        n_syl_in_sequence = [1 1 1 1 1 1];
        
        baseline_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\bl26bl27\bl26bl27_baseline';
        shift_dir = 'H:\6ohda_headphones_all_labeled_songs\6ohda_pitchshift_birds\bl26bl27\bl26bl27_pitch_shift';
        washout_dir = 'F:\DATA\bl26bl27_data\bl26bl27_labeled_songs\bl26bl27_washout';
%         washout_dir = '';
        
    case 'gr193gr194'
        syls = {'a';'s';'d';'f';'q';'w';'e';'j'};
        n_syl_in_sequence = [1 1 1 1 1 1 1 1];
        
        early_baseline = 1;
        switch early_baseline
            case 0 %Day 1-3 baseline
                baseline_dir = 'C:\DATA\gr193gr194_data\gr193gr194_labeled_song\gr193gr194_combined_baseline';
                use_last_n_baseline_days = 1:3;
            case 1 %Day 4-6 baseline
                baseline_dir = 'C:\DATA\gr193gr194_data\gr193gr194_labeled_song\gr193gr194_labeled_baseline';
                use_last_n_baseline_days = 3;
            case 2 %Day 1-6 baseline
                baseline_dir = 'C:\DATA\gr193gr194_data\gr193gr194_labeled_song\gr193gr194_combined_baseline';
                use_last_n_baseline_days = 3;
        end
        shift_dir = 'C:\DATA\gr193gr194_data\gr193gr194_labeled_song\gr193gr194_labeled_neg100shift';
%         washout_dir = '';
        washout_dir = 'C:\DATA\gr193gr194_data\gr193gr194_labeled_song\gr193gr194_labeled_washout';
        
        
    %Here are the 6-OHDA + White noise birds:
    case 'gr93yw43'
        syls = {'b';'c';'e';'f';'h'};
        n_syl_in_sequence = [1 1 1 1 1];
        if strcmp(cond,'prelesion')
            baseline_dir = 'H:\Varun_WN_expts_2014-15_labeled_songs_backup\gr93yw43_prelesion_whitenoise_subset\gr93yw43_baseline_subset';
            shift_dir = '';
            washout_dir = '';
        elseif strcmp(cond,'postlesion')
            baseline_dir = 'H:\gr93yw43_6ohda_lesion_project\gr93yw43_postlesion_white_noise_expt\gr93yw43_labeled_songs\gr93yw43_baseline';
            shift_dir = '';
            washout_dir = '';
        end 
    case 'gr98gy55'
        syls = {'a';'b';'c';'d'};
        n_syl_in_sequence = [1 1 1 1];
        if strcmp(cond,'prelesion')
            baseline_dir = 'H:\Varun_WN_expts_2014-15_labeled_songs_backup\gr98gy55_prelesion_whitenoise_subset\Baseline';
            shift_dir = '';
            washout_dir = '';
        elseif strcmp(cond,'postlesion')
            baseline_dir = 'H:\gr98gy55_6ohda_lesion_project\gr98gy55_postlesion_white_noise_expt\gr98gy55_labeled_songs\gr98gy55_baseline';
            shift_dir = 'H:\gr98gy55_6ohda_lesion_project\gr98gy55_postlesion_white_noise_expt\gr98gy55_labeled_songs\gr98gy55_syl_b_2569_aboveHits';
            washout_dir = '';
        end    
    case 'pk54pu4'
        syls = {'a';'c';'d';'e';'j'};
        n_syl_in_sequence = [1 1 1 1 1];
        if strcmp(cond,'prelesion')
            baseline_dir = 'H:\Varun_WN_expts_2014-15_labeled_songs_backup\pk54pu4_prelesion_whitenoise_subset\Baseline';
            shift_dir = '';
            washout_dir = '';
        elseif strcmp(cond,'postlesion')
            baseline_dir = 'H:\pk54pu4_6ohda_lesion_project\pk54pu4_postlesion_white_noise_expt\pk54pu4_labeled_songs\pk54pu4_baseline';
            shift_dir = '';
            washout_dir = '';
        end
    case 'lb32rd15'
        syls = {'a';'d';'f';'j';'s'};
        n_syl_in_sequence = [1 1 1 1 1];
        if strcmp(cond,'prelesion')
            baseline_dir = 'H:\lb32rd15_6ohda_lesion_project\lb32rd15_prelesion_white_noise_expt\lb32rd15_labeled_songs\lb32rd15_baseline';
            shift_dir = '';
            washout_dir = '';
        elseif strcmp(cond,'postlesion')
            baseline_dir = 'H:\lb32rd15_6ohda_lesion_project\lb32rd15_postlesion_white_noise_expt\lb32rd15_labeled_songs\lb32rd15_baseline';
            shift_dir = '';
            washout_dir = '';
        end
    case 'pk27pu87'
        syls = {'a';'f';'j';'s'};
        n_syl_in_sequence = [1 1 1 1];
        if strcmp(cond,'prelesion')
            baseline_dir = 'H:\pk27pu87_6ohda_lesion_project\pk27pu87_prelesion_white_noise_expt\pk27pu87_labeled_songs\pk27pu87_baseline';
            shift_dir = '';
            washout_dir = '';
        elseif strcmp(cond,'postlesion')
            baseline_dir = 'H:\pk27pu87_6ohda_lesion_project\pk27pu87_postlesion_white_noise_expt\pk27pu87_labeled_songs\pk27pu87_baseline';
            shift_dir = '';
            washout_dir = '';
        end
        
    %Here are Lukas's single syllable shift experiments:
    case 'bl75rd75'
        syls = {'a';'s';'d';'f';'j';'l';'jk';'jkk';'jkkk';'jkkkk';'jkkkkk';'jkkkkkk';'jkkkkkkk';'jkkkkkkkk';'jkkkkkkkkk';'jkkkkkkkkkk'};                   %List of pitch quantified syllables.
        n_syl_in_sequence = [1 1 1 1 1 1 2 3 4 5 6 7 8 9 10 11];
        baseline_dir = 'C:\DATA\Lukas_single_syl_shift_expts\bl75rd75\bl75rd75_baseline';
        shift_dir = 'C:\DATA\Lukas_single_syl_shift_expts\bl75rd75\bl75rd75_100shift';
        washout_dir = '';
    case 'gr44rd54'
        syls = {'a';'b';'c';'f';'g';'i';'ade';'fde';'fdde';'fdde';'e'};                   %List of pitch quantified syllables.
        n_syl_in_sequence = [1 1 1 1 1 1 2 2 2 3 1];
        baseline_dir = 'C:\DATA\Lukas_single_syl_shift_expts\gr44rd54\gr44rd54_baseline';
        shift_dir = 'C:\DATA\Lukas_single_syl_shift_expts\gr44rd54\gr44rd54_neg100shift';
        washout_dir = '';
    case 'or27rd14'
        syls = {'a';'b';'c';'d';'e';'f';'i';'j';'k';'l'};                   %List of pitch quantified syllables.
        n_syl_in_sequence = [1 1 1 1 1 1 1 1 1 1];
        baseline_dir = 'C:\DATA\Lukas_single_syl_shift_expts\or27rd14\or27rd14_baseline';
        shift_dir = 'C:\DATA\Lukas_single_syl_shift_expts\or27rd14\or27rd14_neg100shift';
        washout_dir = '';
    case 'lb4yw24'
        syls = {'a';'b';'c';'d';'e';'f';'r';'u';'v';'w';'x';'y'};                   %List of pitch quantified syllables.
        n_syl_in_sequence = [1 1 1 1 1 1 1 1 1 1 1 1];
        baseline_dir = 'C:\DATA\Lukas_single_syl_shift_expts\lb4yw24\lb4yw24_baseline';
        shift_dir = 'C:\DATA\Lukas_single_syl_shift_expts\lb4yw24\lb4yw24_shift';
        washout_dir = '';
    case 'rd65wh35'
        syls = {'a';'d';'f';'i';'j';'k';'l'};                   %List of pitch quantified syllables.
        n_syl_in_sequence = [1 1 1 1 1 1 1];
        baseline_dir = 'C:\DATA\Lukas_single_syl_shift_expts\rd65wh35\rd65wh35_baseline';
        shift_dir = 'C:\DATA\Lukas_single_syl_shift_expts\rd65wh35\rd65wh35_100shift';
        washout_dir = '';
    case 'pk7r88'
        syls = {'decc';'aabc';'d';'e';'decc';'a';'b';'f';'g'};                   %List of pitch quantified syllables.
        n_syl_in_sequence = [3 4 1 1 4 1 1 1 1];
        baseline_dir = 'C:\DATA\Lukas_single_syl_shift_expts\pk7r88\pk7r88_baseline';
        shift_dir = 'C:\DATA\Lukas_single_syl_shift_expts\pk7r88\pk7r88_shift';
        washout_dir = '';
end
%End of bird-specific settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expt_stage_dir = {baseline_dir; shift_dir; washout_dir}; 
disp(['Baseline song assumed to be in: ' expt_stage_dir{1}]);
disp(['Shifted song assumed to be in: ' expt_stage_dir{2}]);
disp(['Washout song assumed to be in: ' expt_stage_dir{3}]);
num_syls = numel(syls);
num_stages = size(expt_stage_dir,1);
date_labels = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 2: Calculate mean pitch for each syllable, for all baseline days.
cd(baseline_dir);
[dir, num_dirs] = get_dirs();
base_dates = ''; % date from each directory used as a "base" day against which shifts in cents are measured - displayed in the figure title
pitches = cell(num_syls,num_days);
counter = 1;
for d = 1:num_dirs % THIS PRESUMES that last 3 directories returned by ls, are the last 3 baseline or shift days.  Might not be true, check the command line printouts to make sure!
    cd(dir{d}); disp(pwd);
    for s=1:num_syls
        clear weighted_avg;
        weighted_avg = load_pitches(syls{s},n_syl_in_sequence(s));
        if ~isempty(weighted_avg) % if the summary file for syllable s in directory dir{d} existed and contains at least one syllable pitch
            pitches{s,counter} = [pitches{s,counter} weighted_avg(:,1)']; % all pitches for 1 syllable, on 1 day, at 1 experiment stage
        else
            warning(['No pitches found for syllable ' syls{s} ' in directory ' pwd '. Baseline mean pitch for this syllable will necessarily not include this day.']);
        end
    end
    cd ..
    counter = counter + 1;
end 

cd ..
if ~isempty(shift_dir)
    cd(shift_dir);
    [dir, num_dirs] = get_dirs();
    if (num_dirs + 3) > num_days
        for d = num_dirs-2:num_dirs 
            cd(dir{d}); disp(pwd);
            for s=1:num_syls
                clear weighted_avg;
                weighted_avg = load_pitches(syls{s},n_syl_in_sequence(s));
                if ~isempty(weighted_avg) % if the summary file for syllable s in directory dir{d} existed and contains at least one syllable pitch
                    pitches{s,counter} = [pitches{s,counter} weighted_avg(:,1)']; % all pitches for 1 syllable, on 1 day, at 1 experiment stage
                else
                    warning(['No pitches found for syllable ' syls{s} ' in directory ' pwd '. Baseline mean pitch for this syllable will necessarily not include this day.']);
                end
            end
            cd ..
            counter = counter + 1;
        end
    else
       for d = 1:num_dirs 
            cd(dir{d}); disp(pwd);
            for s=1:num_syls
                clear weighted_avg;
                weighted_avg = load_pitches(syls{s},n_syl_in_sequence(s));
                if ~isempty(weighted_avg) % if the summary file for syllable s in directory dir{d} existed and contains at least one syllable pitch
                    pitches{s,counter} = [pitches{s,counter} weighted_avg(:,1)']; % all pitches for 1 syllable, on 1 day, at 1 experiment stage
                else
                    warning(['No pitches found for syllable ' syls{s} ' in directory ' pwd '. Baseline mean pitch for this syllable will necessarily not include this day.']);
                end
            end
            cd ..
            counter = counter + 1;
        end 
    end
    
end

if ~isempty(washout_dir)
    cd ..
    cd(washout_dir);
    [dir, num_dirs] = get_dirs();
    if (counter + num_dirs - 1) <= num_days

        for d = 1:num_dirs 
            cd(dir{d}); disp(pwd);
            for s=1:num_syls
                clear weighted_avg;
                weighted_avg = load_pitches(syls{s},n_syl_in_sequence(s));
                if ~isempty(weighted_avg) % if the summary file for syllable s in directory dir{d} existed and contains at least one syllable pitch
                    pitches{s,counter} = [pitches{s,counter} weighted_avg(:,1)']; % all pitches for 1 syllable, on 1 day, at 1 experiment stage
                else
                    warning(['No pitches found for syllable ' syls{s} ' in directory ' pwd '. Baseline mean pitch for this syllable will necessarily not include this day.']);
                end
            end
            cd ..
            counter = counter + 1;
        end
    else
        warning('Not enough days to include washout data. No washout data included')
    end
end
end

%Utility function
function [sub_dirs, num_dirs] = get_dirs
    dir_data = dir('.');                          %Get the data for the current directory
    dir_index = [dir_data.isdir];                 %Find the index for directories
    sub_dirs = {dir_data(dir_index).name};        %Get a list of the subdirectories
    valid_index = ~ismember(sub_dirs,{'.','..'}); %Find index of subdirectories that are not '.' or '..'   
    sub_dirs = sub_dirs(valid_index);             %Select subdirectories that are not '.' or '..'
    num_dirs = size(sub_dirs,2);                  %Find how many subdirectories there are
end

function weighted_avg = load_pitches(syl,n_syl_in_sequence)
    if length(syl) == 1 %syllable = 'a', or 'b', etc.
        fname_str = ['*_syl_' syl '.mat'];
        fname = ls(fname_str);
    else %syllable is within specific motif 'abcd', or 'aabb', etc.
        fname_str = ['*_syl_' syl(n_syl_in_sequence) '_number_' num2str(n_syl_in_sequence) '_in_sequence_' syl '*'];
        fname = ls(fname_str);
    end
    try
        load(fname,'weighted_avg');
    catch e
       warning(['Summary file ' fname_str ' does not exist. Using empty weighted_avg vector instead.']);
       weighted_avg = [];
    end
    if ~exist('weighted_avg','var') %if not found in summary file
        warning(['weighted_avg NOT FOUND in summary file in directory: ' pwd '. Using empty weighted_avg vector instead.']);
        weighted_avg = [];
    end
end

