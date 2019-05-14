%Script used to extract data from saved svg file of Sober_Brainard_2009 for
%unlesioned no shift birds.

%These values were taken from the svg file directly:
vals1 = [222.4,504.1 235.2,501.8, 248,494.8, 260.8,510.8 273.7,524.7 286.5,507.8 299.3,490.9 312.2,495.1...
325,501, 337.8,501.9, 350.7,503.3, 363.5,499.5, 376.3,516.6 389.2,496.1, 402,484.3];

vals2 = [222.4,504.1, 235.2,512.6, 248,508.1, 260.8,495.4, 273.7,502.3, 286.5,504.3, 299.3,508.3, 312.2,492.1... 
	325,502.2, 337.8,512.8, 350.7,506.3, 363.5,516.2, 376.3,516.7, 389.2,499.8, 402,500.9];

bird1_x = vals1(1:2:end);
bird1_y = vals1(2:2:end);
bird2_x = vals2(1:2:end);
bird2_y = vals2(2:2:end);

%From the line plotting x-axis, I figured out that 222.4 corresponded to
%zero on the x axis:
bird1_x = bird1_x - 222.4;
bird2_x = bird2_x - 222.4;

%Rescaling to the number of days:
bird1_x = bird1_x/bird1_x(end) * 14;
bird2_x = bird2_x/bird2_x(end) * 14;

%For y-axis, the dashed zero line gives us the coordinates for 0 as 504.2:
bird1_y = bird1_y - 504.2;
bird2_y = bird2_y - 504.2;

%Y-axis is positive going downward in SVG files. So we need to invert the
%axis:
bird1_y = -1*bird1_y;
bird2_y = -1*bird2_y;

%The line for the y-axis shows us that the upper limit is 578 which is 73.8
%away from the zero point. This step returns the values in semitones.
bird1_y = bird1_y / 73.8;
bird2_y = bird2_y / 73.8;

%Double-check that the plot now looks exactly the same as the one from
%which we extracted the data:
figure (1)
hold all
plot(bird1_x, bird1_y,'Color',[0 0 0],'LineWidth',2)
plot(bird2_x, bird2_y, 'Color', [0.5 0.5 0.5],'LineWidth',2)
plot(xlim(),[0 0],'k','LineWidth',2)
ylim([-0.6 0.6])
xlim([-0.2 14.2])

%Looks great! So now, let's get the mean and SEM to do stats:
end_means = [bird1_y(end-2:end) bird2_y(end-2:end)];

%Decided to make a plot that gets the mean and SEM using the plots we have
%for the two birds. That plot will be made below:

all_data = [bird1_y; bird2_y];
total_means = mean(all_data);
total_sem = std(all_data)/sqrt(2);

figure (2)
hold all
errorbar(bird1_x,total_means,total_sem,'k','LineWidth',2)
plot(xlim(),[0 0],'k','LineWidth',2)
ylim([-0.6 0.6])
xlim([-0.2 14.2])

