% Written by Varun Saravanan, December 2014.

% The purpose of this function is to simplify the task of defining 
% parameters for each bird when a group analysis to be done on several 
% birds. All parameters are defined here and are passed in the form of a
% structure to the plotting function for simple handling.

% In order to add new birds, simply add a new case for the new bird and
% specify the new parameters in the existing switch case. 

% Input: Bird no integer
% Output: Structure with bird parameters

% Example: [P] = load_bird_params(1)

% Each bird and their corresponding number can be identified in the code
% below. They are also identified in the plotting function 
% make_wn_pitch_shift_and_corr_plots.m

function [ bird_params ] = load_bird_params_2018( bird_no , varargin )

%Basic parameters are defined here:
if nargin > 1
    cond = varargin{1};
end

switch (bird_no)
    case 1  %gr90wh83
        %Enter Params
        bird_params.birdname = 'gr90wh83';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\gr90wh83_prelesion_whitenoise';
        bird_params.shift_direction = 1;            %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'b';'c';'d'};       %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [0 1 2 3];  %Distance of stereotypic syllables from target
        bird_params.target_syl = 1;
        bird_params.same_syls = [2 3 4];            %Position of smae type syllables
        bird_params.diff_syls = [];                 %Position of different type syllables
        bird_params.wn_days = [1 2 3];              %White noise days to consider for learning
        bird_params.motif = 'abcd';
        bird_params.targeted_syl = 'a';
        bird_params.last_wn_day = 3;
                
        
    case 2  %gr93yw43
        %Enter Params
        bird_params.birdname = 'gr93yw43';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\gr93yw43_prelesion_whitenoise';
        bird_params.shift_direction = 1;                    %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'b';'c';'e';'f';'h'};           %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-5 -4 -2 -1 0];    %Distance of stereotypic syllables from target (d is a non-quantifiable syllable)
        bird_params.target_syl = 5;
        bird_params.same_syls = [];                         %Position of same type syllables
        bird_params.diff_syls = [1 2 3 4];                  %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'bcefh';
        bird_params.targeted_syl = 'h';
        bird_params.last_wn_day = 3;

        
    case 3  %gr98gy55
        %Enter Params
        bird_params.birdname = 'gr98gy55';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\gr98gy55_prelesion_whitenoise';
        bird_params.shift_direction = -1;                %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'b';'c';'d'};            %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-1 0 1 2];    %Distance of stereotypic syllables from target
        bird_params.target_syl = 2;
%Use if syl 'a' is a different type syllable:
        bird_params.same_syls = [3 4];                   %Position of same type syllables
        bird_params.diff_syls = [1];                     %Position of different type syllables
%Use these params if agreed to call syl 'a' of this bird a same type
%syllable:
%         bird_params.same_syls = [1 3 4];                   %Position of same type syllables
%         bird_params.diff_syls = [];                     %Position of different type syllables
        bird_params.wn_days = [2 3 4];                   %White noise days to consider for learning
        bird_params.motif = 'abcd';
        bird_params.targeted_syl = 'b';
        bird_params.last_wn_day = 4;
        
    case 4  %pk54pu4
        %Enter Params
        bird_params.birdname = 'pk54pu4';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\pk54pu4_prelesion_whitenoise';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'c';'d';'e';'j'};           %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-2 0 1 2 3];   %Distance of stereotypic syllables from target (b is a non-quantifiable syllable)
        bird_params.target_syl = 2;
        bird_params.same_syls = [];                         %Position of same type syllables
        bird_params.diff_syls = [1 3 4 5];                  %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'acdej';
        bird_params.targeted_syl = 'c';
        bird_params.last_wn_day = 3;
        
    case 5  %or57rd126
        %Enter Params
        bird_params.birdname = 'or57rd126';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\or57rd126_prelesion_whitenoise';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'c';'d'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-2 0 1];       %Distance of stereotypic syllables from target (b is a non-quantifiable syllable)
        bird_params.target_syl = 2;
        bird_params.same_syls = [3];                        %Position of same type syllables
        bird_params.diff_syls = [1];                        %Position of different type syllables
        bird_params.wn_days = [2 3 4];                      %White noise days to consider for learning
        bird_params.motif = 'acd';
        bird_params.targeted_syl = 'c';
        bird_params.last_wn_day = 4;
    
    case 6  %gr93yw43_post 6OHDA lesion
        %Enter Params
        bird_params.birdname = 'gr93yw43';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\Post_lesion_files\gr93yw43_postlesion';
        bird_params.shift_direction = 1;                    %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'b';'c';'e';'f';'h'};           %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-5 -4 -2 -1 0];    %Distance of stereotypic syllables from target (d is a non-quantifiable syllable)
        bird_params.target_syl = 5;
        bird_params.same_syls = [];                         %Position of same type syllables
        bird_params.diff_syls = [1 2 3 4];                  %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'bcefh';
        bird_params.targeted_syl = 'h';
        bird_params.last_wn_day = 3;

        
    case 7  %gr98gy55_post 6OHDA lesion
        %Enter Params
        bird_params.birdname = 'gr98gy55';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\Post_lesion_files\gr98gy55_postlesion';
        bird_params.shift_direction = -1;                %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'b';'c';'d'};            %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-1 0 1 2];    %Distance of stereotypic syllables from target
        bird_params.target_syl = 2;
%Use if syl 'a' is a different type syllable:
        bird_params.same_syls = [3 4];                   %Position of same type syllables
        bird_params.diff_syls = [1];                     %Position of different type syllables
%Use these params if agreed to call syl 'a' of this bird a same type
%syllable:
%         bird_params.same_syls = [1 3 4];                   %Position of same type syllables
%         bird_params.diff_syls = [];                     %Position of different type syllables
        bird_params.wn_days = [1 2 3];                   %White noise days to consider for learning
        bird_params.motif = 'abcd';
        bird_params.targeted_syl = 'b';
        bird_params.last_wn_day = 2;
        
    case 8  %pk54pu4_post 6OHDA lesion
        %Enter Params
        bird_params.birdname = 'pk54pu4';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\Post_lesion_files\pk54pu4_postlesion';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'c';'d';'e'};           %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-2 0 1 2];   %Distance of stereotypic syllables from target
        bird_params.target_syl = 2;
        bird_params.same_syls = [];                         %Position of same type syllables
        bird_params.diff_syls = [1 3 4];                  %Position of different type syllables
        bird_params.wn_days = [4 5 6];                      %White noise days to consider for learning
        bird_params.motif = 'acde';
        bird_params.targeted_syl = 'c';
        bird_params.last_wn_day = 6;
    
    case 9  %lb30rd4
        %Enter Params
        bird_params.birdname = 'lb30rd4';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\lb30rd4_prelesion_whitenoise';
        bird_params.shift_direction = 1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'q';'s';'d';'f'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-3 -3 -2 -1 0];       %Distance of stereotypic syllables from target
        bird_params.target_syl = 5;
        bird_params.same_syls = [3 4];                        %Position of same type syllables
        bird_params.diff_syls = [1 2];                        %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'aqsdf';
        bird_params.targeted_syl = 'f';
        bird_params.last_wn_day = 3;
        
    case 10  %lb32rd15
        %Enter Params
        bird_params.birdname = 'lb32rd15';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\lb32rd15_prelesion_whitenoise';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'j';'j';'a';'d';'d'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-2 -1 0 1 1];       %Distance of stereotypic syllables from target
        bird_params.target_syl = 3;
        bird_params.same_syls = [];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 4 5];                        %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'jjadd';
        bird_params.targeted_syl = 'a';
        bird_params.last_wn_day = 3;
        
    case 11  %pk27pu87 
        %Enter Params
        bird_params.birdname = 'pk27pu87';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\pk27pu87_prelesion_whitenoise';
        bird_params.shift_direction = 1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'s';'s';'f'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-3 -2 -1 0];       %Distance of stereotypic syllables from target
        bird_params.target_syl = 4;
        bird_params.same_syls = [];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 3];                        %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'assf';
        bird_params.targeted_syl = 'f';
        bird_params.last_wn_day = 3;
        
    case 12  %bk26gy38
        %Enter Params
        bird_params.birdname = 'bk26gy38';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\bk26gy38_prelesion_whitenoise';
        bird_params.shift_direction = 1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'s';'d';'f';'j';'k';'l'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-2 -1 0 1 2 3 4];       %Distance of stereotypic syllables from target 
        bird_params.target_syl = 3;
        bird_params.same_syls = [];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 4 5 6 7];                        %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'asdfjkl';
        bird_params.targeted_syl = 'd';
        bird_params.last_wn_day = 3;
    
    case 13  %lb32rd15 - post 6OHDA lesion
        %Enter Params
        bird_params.birdname = 'lb32rd15';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\Post_lesion_files\lb32rd15_postlesion_whitenoise';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'j';'j';'a';'d';'d'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-2 -1 0 1 1];       %Distance of stereotypic syllables from target
        bird_params.target_syl = 3;
        bird_params.same_syls = [];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 4 5];                        %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'jjadd';
        bird_params.targeted_syl = 'a';
        bird_params.last_wn_day = 3;
        
    case 14  %pk27pu87 - post 6OHDA lesion
        %Enter Params
        bird_params.birdname = 'pk27pu87';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\Post_lesion_files\pk27pu87_postlesion_whitenoise';
        bird_params.shift_direction = 1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'s';'s';'f'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-3 -2 -1 0];       %Distance of stereotypic syllables from target
        bird_params.target_syl = 4;
        bird_params.same_syls = [];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 3];                        %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'assf';
        bird_params.targeted_syl = 'f';
        bird_params.last_wn_day = 2;
        
    case 15  %gr90wh83 - post SHAM lesion
        %Enter Params
        bird_params.birdname = 'gr90wh83';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\Post_lesion_files\Sham_lesions\gr90wh83_postlesion_sham';
        bird_params.shift_direction = 1;            %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'b';'c';'d'};       %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [0 1 2 3];  %Distance of stereotypic syllables from target
        bird_params.target_syl = 1;
        bird_params.same_syls = [2 3 4];            %Position of smae type syllables
        bird_params.diff_syls = [];                 %Position of different type syllables
        bird_params.wn_days = [3 4 5];              %White noise days to consider for learning
        bird_params.motif = 'abcd';
        bird_params.targeted_syl = 'a';
        bird_params.last_wn_day = 5;
        
    case 16  %or57rd126 - post SHAM lesion
        %Enter Params
        bird_params.birdname = 'or57rd126';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\Post_lesion_files\Sham_lesions\or57rd126_postlesion_sham';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'c';'d'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-2 0 1];       %Distance of stereotypic syllables from target (b is a non-quantifiable syllable)
        bird_params.target_syl = 2;
        bird_params.same_syls = [3];                        %Position of same type syllables
        bird_params.diff_syls = [1];                        %Position of different type syllables
        bird_params.wn_days = [1 2 3];                      %White noise days to consider for learning
        bird_params.motif = 'acd';
        bird_params.targeted_syl = 'c';
        bird_params.last_wn_day = 3;
        
    case 17  %bk26gy38 - post SHAM lesion
        %Enter Params
        bird_params.birdname = 'bk26gy38';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\Post_lesion_files\Sham_lesions\bk26gy38_postlesion_sham';
        bird_params.shift_direction = 1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'s';'d';'f';'j';'k';'l'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-2 -1 0 1 2 3 4];       %Distance of stereotypic syllables from target 
        bird_params.target_syl = 3;
        bird_params.same_syls = [];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 4 5 6 7];                        %Position of different type syllables
        bird_params.wn_days = [2 3 4];                      %White noise days to consider for learning
        bird_params.motif = 'asdfjkl';
        bird_params.targeted_syl = 'd';
        bird_params.last_wn_day = 4;    
        
    case 18  %lb30rd4 - post SHAM lesion
        %Enter Params
        bird_params.birdname = 'lb30rd4';
        bird_params.top_folder = 'C:\Users\vsarava\Documents\MATLAB\Varun_prelesion_whitenoise_2014\Post_lesion_files\Sham_lesions\lb30rd4_postlesion_sham';
        bird_params.shift_direction = 1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'q';'s';'d';'f'};                   %List of pitch quantified syllables.
        bird_params.syls_away_from_target = [-3 -3 -2 -1 0];       %Distance of stereotypic syllables from target
        bird_params.target_syl = 5;
        bird_params.same_syls = [3 4];                        %Position of same type syllables
        bird_params.diff_syls = [1 2];                        %Position of different type syllables
        bird_params.wn_days = [2 3 4];                      %White noise days to consider for learning
        bird_params.motif = 'aqsdf';
        bird_params.targeted_syl = 'f';
        bird_params.last_wn_day = 4;
       
    %Lukas's headphones generalization birds start here:
    case 19  %pk7r88 - Sam's first single shift experiment
        %Enter Params
        bird_params.birdname = 'pk7r88';
        bird_params.top_folder = 'C:\DATA\Lukas_single_syl_shift_expts\pk7r88';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'decc';'aabc';'d';'e';'decc';'a';'b';'f';'g'};                   %List of pitch quantified syllables.
        bird_params.n_syl_in_sequence = [3 4 1 1 4 1 1 1 1];
        bird_params.syls_away_from_target = [0 -3 -2 -1 1 -5 -4 2 3];       %Distance of stereotypic syllables from target
        bird_params.target_syl = 1;
        bird_params.same_syls = [2 5];                        %Position of same type syllables
        bird_params.diff_syls = [3 4 6 7 8 9];                        %Position of different type syllables
        
    case 20  %lb4yw24 - Lukas's single syl headphones shift experiment
        %Enter Params
        bird_params.birdname = 'lb4yw24';
        bird_params.top_folder = 'C:\DATA\Lukas_single_syl_shift_expts\lb4yw24';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'b';'c';'d';'e';'f';'r';'u';'v';'w';'x';'y'};                   %List of pitch quantified syllables.
        bird_params.n_syl_in_sequence = [1 1 1 1 1 1 1 1 1 1 1 1];
        bird_params.syls_away_from_target = [-1 0 1 2 3 4 5 -2 1 2 3 -3 -2];       %Distance of stereotypic syllables from target
        bird_params.target_syl = 2;
        bird_params.same_syls = [1 3 4 5 6];                        %Position of same type syllables
        bird_params.diff_syls = [7 8 9 10 11 12];
        
    case 21  %or27rd14 - Lukas's single syl headphones shift experiment
        %Enter Params
        bird_params.birdname = 'or27rd14';
        bird_params.top_folder = 'C:\DATA\Lukas_single_syl_shift_expts\or27rd14';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'b';'c';'d';'e';'f';'i';'j';'k';'l'};                   %List of pitch quantified syllables.
        bird_params.n_syl_in_sequence = [1 1 1 1 1 1 1 1 1 1];
        bird_params.syls_away_from_target = [-3 2 -5 -4 -3 3 -2 -1 0 1];       %Distance of stereotypic syllables from target
        bird_params.target_syl = 9;
        bird_params.same_syls = [8 10];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 3 4 5 6 7];
        
    case 22  %gr44rd54 - Lukas's single syl headphones shift experiment
        %Enter Params
        bird_params.birdname = 'gr44rd54';
        bird_params.top_folder = 'C:\DATA\Lukas_single_syl_shift_expts\gr44rd54';
        bird_params.shift_direction = -1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'b';'c';'f';'g';'i';'ade';'fde';'fdde';'fdde';'e'};                   %List of pitch quantified syllables.
        bird_params.n_syl_in_sequence = [1 1 1 1 1 1 2 2 2 3 1];
        bird_params.syls_away_from_target = [-2 2 1 -2 1 1 -1 -1 -1 0 0];       %Distance of stereotypic syllables from target
        bird_params.target_syl = [9 10];
        bird_params.same_syls = [7 8 9];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 3 4 5 6];
        
    case 23  %rd65wh35 - Lukas's single syl headphones shift experiment
        %Enter Params
        bird_params.birdname = 'rd65wh35';
        bird_params.top_folder = 'C:\DATA\Lukas_single_syl_shift_expts\rd65wh35';
        bird_params.shift_direction = 1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'d';'f';'i';'j';'k';'l'};                   %List of pitch quantified syllables.
        bird_params.n_syl_in_sequence = [1 1 1 1 1 1 1];
        bird_params.syls_away_from_target = [1 -2 2 -1 0 0 1];       %Distance of stereotypic syllables from target
        bird_params.target_syl = [5 6];
        bird_params.same_syls = [];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 3 4 7];
        
    case 24  %bl75rd75 - Lukas's single syl headphones shift experiment
        %Enter Params
        bird_params.birdname = 'bl75rd75';
        bird_params.top_folder = 'C:\DATA\Lukas_single_syl_shift_expts\lb75rd75';
        bird_params.shift_direction = 1;                   %+1 for below hits and upward pitch shift. -1 for Above hits and downward pitch shift.
        bird_params.syls = {'a';'s';'d';'f';'j';'l';'jk';'jkk';'jkkk';'jkkkk';'jkkkkk';'jkkkkkk';'jkkkkkkk';'jkkkkkkkk';'jkkkkkkkkk';'jkkkkkkkkkk'};                   %List of pitch quantified syllables.
        bird_params.n_syl_in_sequence = [1 1 1 1 1 1 2 3 4 5 6 7 8 9 10 11];
        bird_params.syls_away_from_target = [-7 -6 -5 -4 -3 1 -2 -1 0 1 2 3 4 5 6 7];       %Distance of stereotypic syllables from target
        bird_params.target_syl = [9];
        bird_params.same_syls = [7 8 10 11 12 13 14 15 16];                        %Position of same type syllables
        bird_params.diff_syls = [1 2 3 4 5 6];
        
    %6-OHDA headphones experiments start here:
    case 25 %bl36yw96 neg100 shift
        bird_params.birdname = 'bl36yw96';
        bird_params.shift_direction = -1;
        bird_params.syls = {'s';'d';'e';'f';'j';'l'};
        
    case 26 %bk9gy9 neg100 shift
        bird_params.birdname = 'bk9gy9';
        bird_params.shift_direction = -1;
        bird_params.syls = {'a';'s';'d';'j';'k';'l'};
        
    case 27 %lb160rd190 neg100 shift
        bird_params.birdname = 'lb160rd190';
        bird_params.shift_direction = -1;
        bird_params.syls = {'s';'d';'f';'e';'h';'j';'k';'t';'l';'i';'a';'g'};
        
    case 28 %bl47yw80 100 shift
        bird_params.birdname = 'bl47yw80';
        bird_params.shift_direction = 1;
        bird_params.syls = {'a';'s';'d';'f';'j';'k';'l';'i'};
        
    case 29 %lb199gy78 100 shift
        bird_params.birdname = 'lb199gy78';
        bird_params.shift_direction = 1;
        bird_params.syls = {'a';'s';'d';'f';'i'};
        
    case 30 %gy46pu6 100 shift
        bird_params.birdname = 'gy46pu6';
        bird_params.shift_direction = 1;
        bird_params.syls = {'a';'s';'f';'g'};
        
    case 31 %bl26bl27 100 shift
        bird_params.birdname = 'bl26bl27';
        bird_params.shift_direction = 1;
        bird_params.syls = {'a';'s';'d';'f';'q';'w'};
        
    case 32 %bl20gr152 no headphones
        bird_params.birdname = 'bl20gr152';
        bird_params.shift_direction = 0;
        bird_params.syls = {'w';'e';'a';'s';'d'};
        
    case 33 %bl23lb23 no headphones
        bird_params.birdname = 'bl23lb23';
        bird_params.shift_direction = 0;
        bird_params.syls = {'a';'s';'d';'f';'r';'w'};
        
    case 34 %pu63pu64 no headphones
        bird_params.birdname = 'pu63pu64';
        bird_params.shift_direction = 0;
        bird_params.syls = {'w';'e';'a';'s';'d'};
        
    case 35 %pu81pu82 no shift
        bird_params.birdname = 'pu81pu82';
        bird_params.shift_direction = 0;
        bird_params.syls = {'a';'s';'d';'f';'g'};
        
    case 36 %bk21bk22 no shift
        bird_params.birdname = 'bk21bk22';
        bird_params.shift_direction = 0;
        bird_params.syls = {'a';'s';'e';'r'};
        
    case 37 %lb140yw4 no shift
        bird_params.birdname = 'lb140yw4';
        bird_params.shift_direction = 0;
        bird_params.syls = {'a';'s';'d';'f';'g'};
        
    case 38 %bk24gy24 no shift
        bird_params.birdname = 'bk24gy24';
        bird_params.shift_direction = 0;
        bird_params.syls = {'a';'s';'d';'k';'q';'w'};
        
    case 39 %lb24or124 no shift
        bird_params.birdname = 'lb24or124';
        bird_params.shift_direction = 0;
        bird_params.syls = {'a';'s';'d';'f';'q'};
        
    case 40 %gr193gr194 neg shift
        bird_params.birdname = 'gr193gr194';
        bird_params.shift_direction = -1;
        bird_params.syls = {'a';'s';'d';'f';'q';'w';'e';'j'};
        
    %Non-lesioned headphones birds start here:
    case 41 %g18g8 neg shift
        bird_params.birdname = 'g18g8_-';
        bird_params.shift_direction = -1;
        bird_params.top_folder = 'C:\Users\vsarava\Documents\Consolidated_Data_for_Python\birddata_100\birddata_100';
        
    case 42 %g44o23 pos shift
        bird_params.birdname = 'g44o23_';
        bird_params.shift_direction = 1;
        bird_params.top_folder = 'C:\Users\vsarava\Documents\Consolidated_Data_for_Python\birddata_100\birddata_100';
        
    case 43 %pk7r88 pos shift
        bird_params.birdname = 'pk7r88_1';
        bird_params.shift_direction = 1;
        bird_params.top_folder = 'C:\Users\vsarava\Documents\Consolidated_Data_for_Python\birddata_100\birddata_100';
        
    case 44 %pk7r88 neg shift
        bird_params.birdname = 'pk7r88_-';
        bird_params.shift_direction = -1;
        bird_params.top_folder = 'C:\Users\vsarava\Documents\Consolidated_Data_for_Python\birddata_100\birddata_100';
        
    case 45 %r12r11 pos shift
        bird_params.birdname = 'r12r11_1';
        bird_params.shift_direction = 1;
        bird_params.top_folder = 'C:\Users\vsarava\Documents\Consolidated_Data_for_Python\birddata_100\birddata_100';
        
    case 46 %r12r11 neg shift
        bird_params.birdname = 'r12r11_-';
        bird_params.shift_direction = -1;
        bird_params.top_folder = 'C:\Users\vsarava\Documents\Consolidated_Data_for_Python\birddata_100\birddata_100';
    
    otherwise
        error('Invalid bird no %d or parameters for bird %d have not been defined',bird_no, bird_no);
end        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following lines of code get the data for each bird and calculate the
%pitch shift of each syllable in cents when compared to the mean pitch over
%all baseline days.


%Adding extra provisions for also including raw data in the form of pitch
%shift in cents wrt baseline mean in the following categories - same type
%and different type as well as distance from target.
    if bird_no < 19
        cd(bird_params.top_folder)
        num_syls = numel(bird_params.motif);
        files = cellstr(ls('*.mat'));
        num_sets = size(files,1);
        bird_params.mean_baseline_pitches = zeros(1,num_syls);
        position = find(bird_params.syls_away_from_target == 0);
        num_days = 6;
        bird_params.data = cell(num_days,num_syls);
        
        Data(1) = load(files{1});
        for j = 1:num_syls
            temp = Data(1).all_weighted_avg(:,j);
            ind = find(temp==0);
            temp(ind) = [];
            bird_params.mean_baseline_pitches(1,j) = mean(temp);
        end
        
        counter = 0;
        for i = 2:num_sets
            Data(i) = load(files{i});
            if sum(bird_params.wn_days==(i-1))==0
                counter = counter + 1;
                continue;
            end
            for j = 1:num_syls
                temp = Data(i).all_weighted_avg(:,j);
                ind = find(temp==0);
                temp(ind) = [];
                bird_params.data_hz{i+2-counter,j} = temp;
                temp = 12*log2(temp/bird_params.mean_baseline_pitches(1,j));
                bird_params.data{i+2-counter,j} = bird_params.shift_direction*temp;
            end
        end
        
        cd Baseline
        files_b = cellstr(ls('*.mat'));
        num_sets_b = size(files_b,1);
        
        for i = 1:num_sets_b
            Data_b(i) = load(files_b{i});
            for j = 1:num_syls
                temp = Data_b(i).all_weighted_avg(:,j);
                ind = find(temp==0);
                temp(ind) = [];
                bird_params.data_hz{i,j} = temp;
                temp = 12*log2(temp/bird_params.mean_baseline_pitches(1,j));
                bird_params.data{i,j} = bird_params.shift_direction*temp;
            end
        end
    
    elseif bird_no < 25
        num_days = 17;
        bird_params.data_hz = get_full_bird_phones_data_vs(bird_params.birdname,num_days)';
        num_syls = length(bird_params.syls);
        bird_params.mean_baseline_pitches = zeros(1,num_syls);
        for j = 1:size(bird_params.data_hz,2)
            temp = [bird_params.data_hz{1,j} bird_params.data_hz{2,j} bird_params.data_hz{3,j}];
            bird_params.mean_baseline_pitches(1,j) = mean(temp);
            for i = 1:size(bird_params.data_hz,1) %Note: The -1 below is because a neg pitch shift causes a pos adaptive change.
                bird_params.data{i,j} = -1*bird_params.shift_direction*12*log2(bird_params.data_hz{i,j}/bird_params.mean_baseline_pitches(1,j));
            end
        end
        
    elseif bird_no < 41
        if exist('cond','var')
            if strcmp(cond,'washout_no_shift')
                num_days = 19;
            else
                num_days = 24;
            end
        else
            num_days = 17;
        end
        bird_params.data_hz = get_full_bird_phones_data_vs(bird_params.birdname,num_days)';
        num_syls = length(bird_params.syls);
        bird_params.mean_baseline_pitches = zeros(1,num_syls);
        for j = 1:size(bird_params.data_hz,2)
            temp = [bird_params.data_hz{1,j} bird_params.data_hz{2,j} bird_params.data_hz{3,j}];
            bird_params.mean_baseline_pitches(1,j) = mean(temp);
            for i = 1:size(bird_params.data_hz,1) %Note: Only normalizing without adjusting for adaptive change.
                bird_params.data{i,j} = 12*log2(bird_params.data_hz{i,j}/bird_params.mean_baseline_pitches(1,j));
            end
        end
        bird_params.data_new = bird_params.data(17:end,:);
        if exist('cond','var')
            if strcmp(cond,'washout')
                bird_params.data = bird_params.data(17:end,:);
                bird_params.mean_end_pitches = zeros(1,num_syls);
                for j = 1:size(bird_params.data,2)
                    bird_params.mean_end_pitches(1,j) = mean(bird_params.data{1,j});
                    for i = 1:size(bird_params.data,1) %Note: Removing mean for each syllable from all days.
                        bird_params.data{i,j} = bird_params.data{i,j} - bird_params.mean_end_pitches(1,j);
                    end
                end
            end
        end
        
    else
        cd(bird_params.top_folder)
        files = cellstr(ls(strcat(bird_params.birdname,'*')));
        temp = load(files{1});
        bird_params.data = [temp.pitch_as_fraction_shift; temp.pitch_as_fraction_washout];
        last_shift_day = size(temp.pitch_as_fraction_shift,1);
        num_syls = size(temp.pitch_as_fraction_shift,2);
        for j = 1:size(bird_params.data,2)
            for i = 1:size(bird_params.data,1) %Note: Converting fractions into semitones
                bird_params.data{i,j} = 12*log2(bird_params.data{i,j});
            end
        end
        if exist('cond','var')
            if ~strcmp(cond,'washout_full')
                bird_params.data = bird_params.data(last_shift_day:(last_shift_day+6),:);
                bird_params.data_unsubtracted = bird_params.data;
                bird_params.mean_end_pitches = zeros(1,num_syls);
                for j = 1:size(bird_params.data,2)
                    bird_params.mean_end_pitches(1,j) = mean(bird_params.data{1,j});
                    for i = 1:size(bird_params.data,1) %Note: Removing mean for each syllable from all days.
                        bird_params.data{i,j} = bird_params.data{i,j} - bird_params.mean_end_pitches(1,j);
                    end
                end
            else
                bird_params.data_new = bird_params.data(last_shift_day:(last_shift_day+6),:);
                bird_params.data = [temp.pitch_as_fraction_shift(1:14,:); temp.pitch_as_fraction_washout];
                for j = 1:size(bird_params.data,2)
                    for i = 1:size(bird_params.data,1) %Note: Converting fractions into semitones
                        bird_params.data{i,j} = 12*log2(bird_params.data{i,j});
                    end
                end
            end
        else
            bird_params.data = bird_params.data(1:14,:);
        end
    end
end