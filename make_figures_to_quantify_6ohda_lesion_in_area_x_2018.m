% Function written June 2015 by Sam Sober to analyze optical densities of
% background, striatum and Area X in TH-stained sections after 6-OHDA or
% sham lesion of Area X. Requires spreadsheets made by the ImageJ macro 
% "quantify_6ohda_lesion_in_area_x" along with added Excel formulas on that
% data.
function make_figures_to_quantify_6ohda_lesion_in_area_x()

current_dir=pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% START HARDCODING SECTION

%This file requires all Excel files with summary data to be stored in the
%same folder. Set the path to that folder in the line below:
cd('C:\2018_lesion_quantification\2018_headphones_birds\6ohda_no_shift_birds\Summary_files')

%Originally, this analysis was performed several different ways and so, we
%retained the structure to do so. However, now you can choose to focus on
%just one "figure_set" and ignore or delete the rest.
figure_set = 4;

if figure_set == 1
    %Load the values from each relevant Excel file:
    [NUM_gr93yw43,TXT_gr93yw43,RAW_gr93yw43] =    xlsread('gr93yw43 6ohda ROI-results');
    [NUM_bk26gy38,TXT_bk26gy38,RAW_bk26gy38] =    xlsread('bk26gy38 sham ROI-results');
    [NUM_or57rd126,TXT_or57rd126,RAW_or57rd126] = xlsread('or57rd126 sham ROI-results');
    [NUM_pk27pu87,TXT_pk27pu87,RAW_pk27pu87] =    xlsread('pk27pu87 6ohda ROI-results');
    [NUM_gr90wh83,TXT_gr90wh83,RAW_gr90wh83]=     xlsread('gr90wh83 sham ROI-results');
    [NUM_pk54pu4,TXT_pk54pu4,RAW_pk54pu4] =       xlsread('pk54pu4 6ohda ROI-results');
    [NUM_gr98gy55,TXT_gr98gy55,RAW_gr98gy55] =    xlsread('gr98gy55 6ohda ROI-results');
    [NUM_lb32rd15,TXT_lb32rd15,RAW_lb32rd15] =    xlsread('lb32rd15 6ohda ROI-results');
    [NUM_lb30rd4,TXT_lb30rd4,RAW_lb30rd4] =       xlsread('lb30rd4 sham ROI-results');
    [NUM_bk21bk22,TXT_bk21bk22,RAW_bk21bk22] =    xlsread('bk21bk22 6ohda ROI-results');
    [NUM_br24gy24,TXT_br24gy24,RAW_br24gy24] =    xlsread('br24gy24 6ohda ROI-results');
    [NUM_lb24or124,TXT_lb24or124,RAW_lb24or124] = xlsread('lb24or124 6ohda ROI-results');
    [NUM_lb140yw4,TXT_lb140yw4,RAW_lb140yw4] =    xlsread('lb140yw4 6ohda ROI-results');
    [NUM_pu81pu82,TXT_pu81pu82,RAW_pu81pu82] =    xlsread('pu81pu82 6ohda ROI-results');
    titlestr = 'unnormalized';
elseif figure_set == 2
    [NUM_pu81pu82,TXT_pu81pu82,RAW_pu81pu82] =    xlsread('gr93yw43 6ohda normalized background ROI_results');
    [NUM_bk26gy38,TXT_bk26gy38,RAW_bk26gy38] =    xlsread('bk26gy38 sham normalized background ROI-results');
    [NUM_or57rd126,TXT_or57rd126,RAW_or57rd126] = xlsread('or57rd126 sham normalized background ROI-results');
    [NUM_bk21bk22,TXT_bk21bk22,RAW_bk21bk22] =    xlsread('pk27pu87 6ohda normalized background ROI-results');
    [NUM_gr90wh83,TXT_gr90wh83,RAW_gr90wh83]=     xlsread('gr90wh83 sham normalized background ROI-results');
    [NUM_br24gy24,TXT_br24gy24,RAW_br24gy24] =       xlsread('pk54pu4 6ohda normalized background ROI-results');
    [NUM_lb24or124,TXT_lb24or124,RAW_lb24or124] =    xlsread('gr98gy55 6ohda normalized background ROI-results');
    [NUM_lb140yw4,TXT_lb140yw4,RAW_lb140yw4] =    xlsread('lb32rd15 6ohda normalized background ROI-results');
    [NUM_lb30rd4,TXT_lb30rd4,RAW_lb30rd4] =       xlsread('lb30rd4 sham normalized background ROI-results'); 
    titlestr = 'background-normalized';
elseif figure_set == 3
    [NUM_pu81pu82,TXT_pu81pu82,RAW_pu81pu82] =    xlsread('gr93yw43 6ohda normalized striatum ROI_results');
    [NUM_bk26gy38,TXT_bk26gy38,RAW_bk26gy38] =    xlsread('bk26gy38 sham normalized striatum ROI-results');
    [NUM_or57rd126,TXT_or57rd126,RAW_or57rd126] = xlsread('or57rd126 sham normalized striatum ROI-results');
    [NUM_bk21bk22,TXT_bk21bk22,RAW_bk21bk22] =    xlsread('pk27pu87 6ohda normalized striatum ROI-results');
    [NUM_gr90wh83,TXT_gr90wh83,RAW_gr90wh83]=     xlsread('gr90wh83 sham normalized striatum ROI-results');
    [NUM_br24gy24,TXT_br24gy24,RAW_br24gy24] =       xlsread('pk54pu4 6ohda normalized striatum ROI-results');
    [NUM_lb24or124,TXT_lb24or124,RAW_lb24or124] =    xlsread('gr98gy55 6ohda normalized striatum ROI-results');
    [NUM_lb140yw4,TXT_lb140yw4,RAW_lb140yw4] =    xlsread('lb32rd15 6ohda normalized striatum ROI-results');
    [NUM_lb30rd4,TXT_lb30rd4,RAW_lb30rd4] =       xlsread('lb30rd4 sham normalized striatum ROI-results');
    titlestr = 'striatum-normalized';
elseif figure_set == 4
    [NUM_gr93yw43,TXT_gr93yw43,RAW_gr93yw43] =    xlsread('gr93yw43 6ohda ROI-results');
    [NUM_bk26gy38,TXT_bk26gy38,RAW_bk26gy38] =    xlsread('bk26gy38 sham ROI-results');
    [NUM_or57rd126,TXT_or57rd126,RAW_or57rd126] = xlsread('or57rd126 sham ROI-results');
    [NUM_pk27pu87,TXT_pk27pu87,RAW_pk27pu87] =    xlsread('pk27pu87 6ohda ROI-results');
    [NUM_gr90wh83,TXT_gr90wh83,RAW_gr90wh83]=     xlsread('gr90wh83 sham ROI-results');
    [NUM_pk54pu4,TXT_pk54pu4,RAW_pk54pu4] =       xlsread('pk54pu4 6ohda ROI-results');
    [NUM_gr98gy55,TXT_gr98gy55,RAW_gr98gy55] =    xlsread('gr98gy55 6ohda ROI-results');
    [NUM_lb32rd15,TXT_lb32rd15,RAW_lb32rd15] =    xlsread('lb32rd15 6ohda ROI-results');
    [NUM_lb30rd4,TXT_lb30rd4,RAW_lb30rd4] =       xlsread('lb30rd4 sham ROI-results');
    [NUM_bk9gy9,TXT_bk9gy9,RAW_bk9gy9] =    xlsread('bk9gy9 6ohda ROI-results');
    [NUM_bl26bl27,TXT_bl26bl27,RAW_bl26bl27] =    xlsread('bl26bl27 6ohda ROI-results');
    [NUM_bl36yw96,TXT_bl36yw96,RAW_bl36yw96] = xlsread('bl36yw96 6ohda ROI-results');
    [NUM_bl47yw80,TXT_bl47yw80,RAW_bl47yw80] =    xlsread('bl47yw80 6ohda ROI-results');
    [NUM_gr193gr194,TXT_gr193gr194,RAW_gr193gr194] =    xlsread('gr193gr194 6ohda ROI-results');
    [NUM_gy46pu6,TXT_gy46pu6,RAW_gy46pu6] =    xlsread('gy46pu6 6ohda ROI-results');
    [NUM_lb160rd190,TXT_lb160rd190,RAW_lb160rd190] =    xlsread('lb160rd190 6ohda ROI-results');
    [NUM_lb199gy78,TXT_lb199gy78,RAW_lb199gy78] =    xlsread('lb199gy78 6ohda ROI-results');
    [NUM_bk21bk22,TXT_bk21bk22,RAW_bk21bk22] =    xlsread('bk21bk22 6ohda ROI-results');
    [NUM_br24gy24,TXT_br24gy24,RAW_br24gy24] =    xlsread('br24gy24 6ohda ROI-results');
    [NUM_lb24or124,TXT_lb24or124,RAW_lb24or124] = xlsread('lb24or124 6ohda ROI-results');
    [NUM_lb140yw4,TXT_lb140yw4,RAW_lb140yw4] =    xlsread('lb140yw4 6ohda ROI-results');
    [NUM_pu81pu82,TXT_pu81pu82,RAW_pu81pu82] =    xlsread('pu81pu82 6ohda ROI-results');
    [NUM_bl20gr152,TXT_bl20gr152,RAW_bl20gr152] = xlsread('bl20gr152 6ohda ROI-results');
    [NUM_bl23lb23,TXT_bl23lb23,RAW_bl23lb23] =    xlsread('bl23lb23 6ohda ROI-results');
    [NUM_pu63pu64,TXT_pu63pu64,RAW_pu63pu64] =    xlsread('pu63pu64 6ohda ROI-results');
    titlestr = 'unnormalized';
end

%For each bird used, set the bird name, a color for plotting and a symbol
%chosen for the individual plots: (Note, if only using the cdf plot, one
%can avoid setting the colors and symbols.
birdname_arr{1}='gr93yw43';col_vec(1)='r';symb_vec(1)='.';    % 6OHDA
birdname_arr{2}='bk26gy38';col_vec(2)='k';symb_vec(2)='o';    % SHAM
birdname_arr{3}='or57rd126';col_vec(3)='k';symb_vec(3)='s';   % SHAM
birdname_arr{4}='pk27pu87';col_vec(4)='g';symb_vec(4)='.';    % 6OHDA 
birdname_arr{5}='gr90wh83';col_vec(5)='k';symb_vec(5)='d';    % SHAM
birdname_arr{6}='pk54pu4';col_vec(6)='b';symb_vec(6)='.';     % 6OHDA 
birdname_arr{7}='gr98gy55';col_vec(7)='m';symb_vec(7)='.';    % 6OHDA 
birdname_arr{8}='lb32rd15';col_vec(8)='c';symb_vec(8)='.';    % 6OHDA 
birdname_arr{9}='lb30rd4';col_vec(9)='k';symb_vec(9)='x';    % SHAM

if figure_set == 1
    birdname_arr{10}='bk21bk22';col_vec(10)='r';symb_vec(10)='.';    % 6OHDA
    birdname_arr{11}='br24gy24';col_vec(11)='g';symb_vec(11)='.';     % 6OHDA 
    birdname_arr{12}='lb24or124';col_vec(12)='b';symb_vec(12)='.';    % 6OHDA 
    birdname_arr{13}='lb140yw4';col_vec(13)='m';symb_vec(13)='.';    % 6OHDA 
    birdname_arr{14}='pu81pu82';col_vec(14)='c';symb_vec(14)='.';    % 6OHDA
    ohda_vector=[10 11 12 13 14];
elseif figure_set == 4
    birdname_arr{10}='bk9gy9';col_vec(10)='r';symb_vec(10)='.';    % 6OHDA
    birdname_arr{11}='bl26bl27';col_vec(11)='g';symb_vec(11)='.';     % 6OHDA 
    birdname_arr{12}='bl36yw96';col_vec(12)='b';symb_vec(12)='.';    % 6OHDA 
    birdname_arr{13}='bl47yw80';col_vec(13)='m';symb_vec(13)='.';    % 6OHDA 
    birdname_arr{14}='gr193gr194';col_vec(14)='c';symb_vec(14)='.';    % 6OHDA
    birdname_arr{15}='gy46pu6';col_vec(15)='g';symb_vec(15)='x';    % 6OHDA 
    birdname_arr{16}='lb160rd190';col_vec(16)='r';symb_vec(16)='x';    % 6OHDA 
    birdname_arr{17}='lb199gy78';col_vec(17)='b';symb_vec(17)='x';    % 6OHDA
    birdname_arr{18}='bk21bk22';col_vec(10)='r';symb_vec(10)='.';    % 6OHDA
    birdname_arr{19}='br24gy24';col_vec(11)='g';symb_vec(11)='.';     % 6OHDA 
    birdname_arr{20}='lb24or124';col_vec(12)='b';symb_vec(12)='.';    % 6OHDA 
    birdname_arr{21}='lb140yw4';col_vec(13)='m';symb_vec(13)='.';    % 6OHDA 
    birdname_arr{22}='pu81pu82';col_vec(14)='c';symb_vec(14)='.';    % 6OHDA
    birdname_arr{23}='bl20gr152';col_vec(12)='b';symb_vec(12)='.';    % 6OHDA 
    birdname_arr{24}='bl23lb23';col_vec(13)='m';symb_vec(13)='.';    % 6OHDA 
    birdname_arr{25}='pu63pu64';col_vec(14)='c';symb_vec(14)='.';    % 6OHDA
    ohda_vector=[10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
end

%For plotting individual sham birds in varying shades of grayscale:
sham_grey_col_array{2}=[0 0 0];
sham_grey_col_array{3}=[.3 .3 .3];
sham_grey_col_array{5}=[.6 .6 .6];
sham_grey_col_array{9}=[.9 .9 .9];

% IDs of shams and lesions
sham_vector=[2 3 5 9];
old_ohda_vector=[1 4 6 7 8];


% NOTE - all row numberings need to be checked here, and checked separately for each bird
for x=1:length(birdname_arr)
    %Load the relevant variables into new vectors that will be used for
    %analysis:
    eval(sprintf('%s_raw_ratio_vec=NUM_%s(25,:);',birdname_arr{x},birdname_arr{x}));
    eval(sprintf('%s_BL_subtracted_ratio_vec=NUM_%s(26,:);',birdname_arr{x},birdname_arr{x}));
    eval(sprintf('%s_X_minux_BG_vec=NUM_%s(17,:);',birdname_arr{x},birdname_arr{x}));
    eval(sprintf('%s_Str_minux_BG_vec=NUM_%s(14,:);',birdname_arr{x},birdname_arr{x}));

    eval(sprintf('bg_raw_%s=NUM_%s(1,:);',birdname_arr{x},birdname_arr{x}));
    eval(sprintf('X_raw_%s=NUM_%s(16,:);',birdname_arr{x},birdname_arr{x}));
    eval(sprintf('Str_raw_%s=NUM_%s(13,:);',birdname_arr{x},birdname_arr{x}));
end

% END HARDCODING SECTION
%%%%%%%%%%%%%%%%%%%%%%%%

raw_ratio_vec_COMBINED_6OHDA=[];
raw_ratio_vec_COMBINED_old_6OHDA=[];
raw_ratio_vec_COMBINED_SHAM=[];
%raw_X_vec_COMBINED_SHAM=[];

%Create a vector that has the combined raw ratio for each subset of birds
%which will be used to produce the cdf plot. Note a potential downside of
%this analysis is that every section is given equal weight regardless of
%the size of Area X within the section:
for x=1:length(sham_vector)
    eval(sprintf('raw_ratio_vec_COMBINED_SHAM=[raw_ratio_vec_COMBINED_SHAM %s_raw_ratio_vec];',birdname_arr{sham_vector(x)}));
    %eval(sprintf('raw_X_vec_COMBINED_SHAM=[raw_X_vec_COMBINED_SHAM X_raw_%s]',birdname_arr{sham_vector(x)}));
end
for x=1:length(ohda_vector)
    eval(sprintf('raw_ratio_vec_COMBINED_6OHDA=[raw_ratio_vec_COMBINED_6OHDA %s_raw_ratio_vec];',birdname_arr{ohda_vector(x)}));
end
for x=1:length(old_ohda_vector)
    eval(sprintf('raw_ratio_vec_COMBINED_old_6OHDA=[raw_ratio_vec_COMBINED_old_6OHDA %s_raw_ratio_vec];',birdname_arr{old_ohda_vector(x)}));
end

%Generate the plot with the cdfs for each group:
figure('units','normalized','outerposition',[0 0 1 1]); % full screen
hold on
a=cdfplot(raw_ratio_vec_COMBINED_SHAM);set(a,'color','k','linew',2)
b=cdfplot(raw_ratio_vec_COMBINED_6OHDA);set(b,'color','r','linew',2)
%c=cdfplot(raw_ratio_vec_COMBINED_old_6OHDA);set(c,'color','b','linew',2)
% title(['Black: sham (n=' num2str(numel(sham_vector)) ') Blue: Old 6OHDA (n= ' num2str(numel(old_ohda_vector)) '), Red: 6OHDA (n= ' num2str(numel(ohda_vector)) '), all figs ' titlestr]);xlabel('Ratio of Raw\_X/Raw\_Str')
title(['Black: sham (n=' num2str(numel(sham_vector)) ')  Red: 6OHDA (n= ' num2str(numel(ohda_vector)) '), all figs ' titlestr]);xlabel('Ratio of Raw\_X/Raw\_Str')
set(gca,'xtick',[0:.5:2],'ytick',[0:.25:1],'xlim',[0 2],'ylim',[0 1])
%Perform a ks-test to check for differences between the sham group and
%lesioned group:
[H_combo,P_combo,K_combo] =kstest2(raw_ratio_vec_COMBINED_SHAM,raw_ratio_vec_COMBINED_6OHDA)

msize=3;
text_move_frac=0.1;

%Get example points in the cdf plot (note that this has been hardcoded
%based on the example sections chosen. Currently hardcoded for Varun's
%6-OHDA + headphones paper:
s1=scatter(0.9277,0.1556,160,'r','s')
s1.LineWidth = 2;
s1.MarkerEdgeColor = 'r';
s1.MarkerFaceColor = [1 1 1];

s2=scatter(1.1532,0.4951,160,'k')
s2.LineWidth = 2;
s2.MarkerEdgeColor = 'k';
s2.MarkerFaceColor = [1 1 1];

%Shade the area that is beyond the 5th percentile in sham lesioned birds as
%lesioned section. Again, hardcoded for saline birds from Lukas's 2016
%paper.
a = area([1.0467 2], [1 1],'LineStyle',':')
plot([0 1.0467],[0.3746 0.3746],'k','LineWidth',2)

%Varun typically stops execution here. The rest of the figures are for
%various other ways of analysis performed for Lukas's 2016 paper.

figure('units','normalized','outerposition',[0 0 1 1]); % full screen
%Reproduce the same plot as above:
subplot(3,3,1);hold on
a=cdfplot(raw_ratio_vec_COMBINED_SHAM);set(a,'color','k','linew',2)
b=cdfplot(raw_ratio_vec_COMBINED_6OHDA);set(b,'color','r','linew',2)
title(['Black: sham (n=' num2str(numel(sham_vector)) ') Red: 6OHDA (n= ' num2str(numel(ohda_vector)) '), all figs ' titlestr]);xlabel('Ratio of Raw\_X/Raw\_Str')
set(gca,'xtick',[0:.5:2],'ytick',[0:.25:1],'xlim',[0 2],'ylim',[0 1])
[H_combo,P_combo,K_combo] =kstest2(raw_ratio_vec_COMBINED_SHAM,raw_ratio_vec_COMBINED_6OHDA)

%Get individual traces for each lesioned bird with a combined sham trace.
subplot(3,3,2);hold on
a=cdfplot(raw_ratio_vec_COMBINED_SHAM);set(a,'color','k','linew',2)
for x=1:length(ohda_vector)
    eval(sprintf('tmp=cdfplot(%s_raw_ratio_vec);',birdname_arr{ohda_vector(x)}))
    set(tmp,'color',col_vec(ohda_vector(x)),'marker',symb_vec(ohda_vector(x)),'markersize',msize,'linew',2)
end
title('Black: sham; Color: 6OHDA');xlabel('Ratio of Raw\_X/Raw\_Str')
xl=get(gca,'xlim');yl=get(gca,'ylim');
for x=1:length(ohda_vector)
    text(xl(1)+text_move_frac*diff(xl),yl(2)-text_move_frac*x*diff(yl),birdname_arr(ohda_vector(x)),'color',col_vec(ohda_vector(x)),'fontweight','bold','fontsize',12)
    eval(sprintf('[H,P,K] =kstest2(raw_ratio_vec_COMBINED_SHAM,%s_raw_ratio_vec)',birdname_arr{ohda_vector(x)}))
end
set(gca,'xtick',[0:.5:2],'ytick',[0:.25:1],'xlim',[0 2],'ylim',[0 1])

%Get individual traces for both saline and lesioned birds.
subplot(3,3,3);hold on
for x=1:length(sham_vector)
    eval(sprintf('tmp=cdfplot(%s_raw_ratio_vec);',birdname_arr{sham_vector(x)}))
    set(tmp,'color',sham_grey_col_array{sham_vector(x)},'marker',symb_vec(sham_vector(x)),'markersize',msize,'linew',2)
end
for x=1:length(ohda_vector)
    eval(sprintf('tmp=cdfplot(%s_raw_ratio_vec);',birdname_arr{ohda_vector(x)}))
    set(tmp,'color',col_vec(ohda_vector(x)),'marker',symb_vec(ohda_vector(x)),'markersize',msize,'linew',2)
end
title('Black: sham; Color: 6OHDA');xlabel('Ratio of Raw\_X/Raw\_Str')
xl=get(gca,'xlim');yl=get(gca,'ylim');
for x=1:length(ohda_vector)
    text(xl(1)+text_move_frac*diff(xl),yl(2)-text_move_frac*x*diff(yl),birdname_arr(ohda_vector(x)),'color',col_vec(ohda_vector(x)),'fontweight','bold','fontsize',12)
end

%From here on, it is the same analysis but with slightly different metrics
%in each case:
subplot(3,3,4);hold on
for x=1:length(sham_vector)
    eval(sprintf('tmp=cdfplot(%s_raw_ratio_vec/max(%s_raw_ratio_vec));',birdname_arr{sham_vector(x)},birdname_arr{sham_vector(x)}))
    set(tmp,'color',sham_grey_col_array{sham_vector(x)},'marker',symb_vec(sham_vector(x)),'markersize',msize,'linew',2)
end
for x=1:length(ohda_vector)
    eval(sprintf('tmp=cdfplot(%s_raw_ratio_vec/max(%s_raw_ratio_vec));',birdname_arr{ohda_vector(x)},birdname_arr{ohda_vector(x)}))
    set(tmp,'color',col_vec(ohda_vector(x)),'marker',symb_vec(ohda_vector(x)),'markersize',msize,'linew',2)
end
title('Black: sham; Color: 6OHDA');xlabel('Ratio of Raw\_X/Raw\_Str, normalized by GREATEST value of this ratio')
xl=get(gca,'xlim');yl=get(gca,'ylim');
for x=1:length(ohda_vector)
    text(xl(1)+text_move_frac*diff(xl),yl(2)-text_move_frac*x*diff(yl),birdname_arr(ohda_vector(x)),'color',col_vec(ohda_vector(x)),'fontweight','bold','fontsize',12)
end

subplot(3,3,5);hold on
for x=1:length(sham_vector)
    eval(sprintf('tmp=cdfplot(%s_X_minux_BG_vec/max(%s_X_minux_BG_vec));',birdname_arr{sham_vector(x)},birdname_arr{sham_vector(x)}))
    set(tmp,'color',sham_grey_col_array{sham_vector(x)},'marker',symb_vec(sham_vector(x)),'markersize',msize,'linew',2)
end
for x=1:length(ohda_vector)
    eval(sprintf('tmp=cdfplot(%s_X_minux_BG_vec/max(%s_X_minux_BG_vec));',birdname_arr{ohda_vector(x)},birdname_arr{ohda_vector(x)}))
    set(tmp,'color',col_vec(ohda_vector(x)),'marker',symb_vec(ohda_vector(x)),'markersize',msize,'linew',2)
end
title('Black: sham; Color: 6OHDA');xlabel('BG-subtracted\_X, normalized by GREATEST value')
xl=get(gca,'xlim');yl=get(gca,'ylim');
for x=1:length(ohda_vector)
    text(xl(1)+text_move_frac*diff(xl),yl(2)-text_move_frac*x*diff(yl),birdname_arr(ohda_vector(x)),'color',col_vec(ohda_vector(x)),'fontweight','bold','fontsize',12)
end

subplot(3,3,6);hold on
for x=1:length(sham_vector)
    eval(sprintf('tmp=cdfplot(%s_Str_minux_BG_vec/max(%s_Str_minux_BG_vec));',birdname_arr{sham_vector(x)},birdname_arr{sham_vector(x)}))
    set(tmp,'color',sham_grey_col_array{sham_vector(x)},'marker',symb_vec(sham_vector(x)),'markersize',msize,'linew',2)
end
for x=1:length(ohda_vector)
    eval(sprintf('tmp=cdfplot(%s_Str_minux_BG_vec/max(%s_Str_minux_BG_vec));',birdname_arr{ohda_vector(x)},birdname_arr{ohda_vector(x)}))
    set(tmp,'color',col_vec(ohda_vector(x)),'marker',symb_vec(ohda_vector(x)),'markersize',msize,'linew',2)
end
title('Black: sham; Color: 6OHDA');xlabel('BG-subtracted\_Striatum, normalized by GREATEST value')
xl=get(gca,'xlim');yl=get(gca,'ylim');
for x=1:length(ohda_vector)
    text(xl(1)+text_move_frac*diff(xl),yl(2)-text_move_frac*x*diff(yl),birdname_arr(ohda_vector(x)),'color',col_vec(ohda_vector(x)),'fontweight','bold','fontsize',12)
end

% a bit more analysis
raw_X_over_mean_Str_vec_COMBINED_SHAM=[];
raw_X_over_max_X_within_bird_vec_COMBINED_SHAM=[];
raw_Str_over_max_Str_within_bird_vec_COMBINED_SHAM=[];
for x=1:length(sham_vector)
    eval(sprintf('raw_X_over_mean_Str_vec_COMBINED_SHAM=[raw_X_over_mean_Str_vec_COMBINED_SHAM X_raw_%s/mean(Str_raw_%s)];',birdname_arr{sham_vector(x)},birdname_arr{sham_vector(x)}));
    eval(sprintf('raw_X_over_max_X_within_bird_vec_COMBINED_SHAM=[raw_X_over_max_X_within_bird_vec_COMBINED_SHAM X_raw_%s/max(X_raw_%s)];',birdname_arr{sham_vector(x)},birdname_arr{sham_vector(x)}));
    eval(sprintf('raw_Str_over_max_Str_within_bird_vec_COMBINED_SHAM=[raw_Str_over_max_Str_within_bird_vec_COMBINED_SHAM Str_raw_%s/max(Str_raw_%s)];',birdname_arr{sham_vector(x)},birdname_arr{sham_vector(x)}));
end

subplot(3,3,7);hold on
a=cdfplot(raw_X_over_mean_Str_vec_COMBINED_SHAM);set(a,'color','k','linew',2)
for x=1:length(ohda_vector)
    eval(sprintf('tmp=cdfplot(X_raw_%s/mean(Str_raw_%s));',birdname_arr{ohda_vector(x)},birdname_arr{ohda_vector(x)}))
    set(tmp,'color',col_vec(ohda_vector(x)),'marker',symb_vec(ohda_vector(x)),'markersize',msize,'linew',2)
end
title('Black: sham; Color: 6OHDA');xlabel('Ratio of Raw\_X/mean(Raw\_Str)')

subplot(3,3,8);hold on
a=cdfplot(raw_X_over_max_X_within_bird_vec_COMBINED_SHAM);set(a,'color','k','linew',2)
for x=1:length(ohda_vector)
    eval(sprintf('tmp=cdfplot(X_raw_%s/max(X_raw_%s));',birdname_arr{ohda_vector(x)},birdname_arr{ohda_vector(x)}))
    set(tmp,'color',col_vec(ohda_vector(x)),'markersize',msize,'linew',2)
    eval(sprintf('[H,P,K] =kstest2(raw_X_over_max_X_within_bird_vec_COMBINED_SHAM,X_raw_%s/max(X_raw_%s))',birdname_arr{ohda_vector(x)},birdname_arr{ohda_vector(x)}))
end
title('Black: sham; Color: 6OHDA');xlabel('Ratio of Raw\_X/max(Raw\_X)')
set(gca,'xtick',[0:.25:1],'ytick',[0:.25:1],'xlim',[0 1],'ylim',[0 1])

subplot(3,3,9);hold on
a=cdfplot(raw_Str_over_max_Str_within_bird_vec_COMBINED_SHAM);set(a,'color','k','linew',2)
for x=1:length(ohda_vector)
    eval(sprintf('tmp=cdfplot(Str_raw_%s/max(Str_raw_%s));',birdname_arr{ohda_vector(x)},birdname_arr{ohda_vector(x)}))
    set(tmp,'color',col_vec(ohda_vector(x)),'markersize',msize,'linew',2)
end
title('Black: sham; Color: 6OHDA');xlabel('Ratio of Raw\_Str/max(Raw\_Str)')
set(gca,'xtick',[0:.25:1],'ytick',[0:.25:1],'xlim',[0 1],'ylim',[0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('units','normalized','outerposition',[0 0 1 1]); % full screen
w(1)=subplot(3,3,1);hold on;
for x=1:length(birdname_arr)
    eval(sprintf('tmp=cdfplot(%s_raw_ratio_vec);',birdname_arr{x}))
    set(tmp,'color',col_vec(x),'marker',symb_vec(x),'markersize',msize)
end
title(['Ratio of Raw\_X/Raw\_Str, all figs ' titlestr])

w(2)=subplot(3,3,4);hold on;
for x=1:length(birdname_arr)
    eval(sprintf('tmp=cdfplot(%s_BL_subtracted_ratio_vec);',birdname_arr{x}))
    set(tmp,'color',col_vec(x),'marker',symb_vec(x),'markersize',msize)
end
title('Ratio of BL\_subtracted\_X/BL\_subtracted\_Str')
linkaxes(w,'xy');set(gca,'xlim',[0 2])

v(1)=subplot(3,3,2);hold on;
for x=1:length(birdname_arr)
    eval(sprintf('tmp=cdfplot(%s_X_minux_BG_vec);',birdname_arr{x}))
    set(tmp,'color',col_vec(x),'marker',symb_vec(x),'markersize',msize)
end
title('X minus BG')

v(2)=subplot(3,3,5);hold on;
for x=1:length(birdname_arr)
    eval(sprintf('tmp=cdfplot(%s_Str_minux_BG_vec);',birdname_arr{x}))
    set(tmp,'color',col_vec(x),'marker',symb_vec(x),'markersize',msize)
end
title('Str minus BG')
linkaxes(v,'xy');set(gca,'xlim',[0 140])

z(1)=subplot(3,3,3);hold on;
for x=1:length(birdname_arr)
    eval(sprintf('tmp=cdfplot(bg_raw_%s);',birdname_arr{x}))
    set(tmp,'color',col_vec(x),'marker',symb_vec(x),'markersize',msize)
end
title('Raw BG')

z(2)=subplot(3,3,6);
hold on;
for x=1:length(birdname_arr)
    eval(sprintf('tmp=cdfplot(X_raw_%s);',birdname_arr{x}))
    set(tmp,'color',col_vec(x),'marker',symb_vec(x),'markersize',msize)
end
title('Raw X')

z(3)=subplot(3,3,9);
hold on;
for x=1:length(birdname_arr)
    eval(sprintf('tmp=cdfplot(Str_raw_%s);',birdname_arr{x}))
    set(tmp,'color',col_vec(x),'marker',symb_vec(x),'markersize',msize)
end
title('Raw Str')
linkaxes(z,'xy');set(gca,'xlim',[0 160])

z(3)=subplot(3,3,7);hold on;
for x=1:length(birdname_arr)
    eval(sprintf('tmp=cdfplot(X_raw_%s/mean(bg_raw_%s));',birdname_arr{x},birdname_arr{x}))
    set(tmp,'color',col_vec(x),'marker',symb_vec(x),'markersize',msize)
end
title('Raw X/Mean BG')

z(3)=subplot(3,3,8);hold on;
for x=1:length(birdname_arr)
    eval(sprintf('tmp=cdfplot(X_raw_%s/mean(Str_raw_%s));',birdname_arr{x},birdname_arr{x}))
    set(tmp,'color',col_vec(x),'marker',symb_vec(x),'markersize',msize)
end
title('Raw X/Mean Str')

cd(current_dir)
end
