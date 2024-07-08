%% pupiltrack
% This code is used to extract and plot pupil data extracted from
% deeplabcut (DLC) using Justin's 8 point trainer.
% Original author Justin Aronson 10/19
% last update by Angela Nietz 6/13/2022

%%
clear
day_in = uigetdir();  % lets you pick a folder to analyze
dir_in = [day_in,'\eye\']; % tells program to look in eye folder
cd(dir_in); % makes the current directory the folder you picked
mouseID = day_in(31:end); %mouseID
date = day_in(24:29); %experiment date
FileList = dir(fullfile([dir_in,'*filtered.csv'])); %loads point tracking data from DLC
aq_freq = 20 % video aquisition frequency (make sure to change for 20 Hz videos)
tr_length = 330 %trial length
for File = 1:length(FileList); % list of analyzed videos in chosen folder
    raw_data(:,:,File) = importdata(FileList(File).name); % imports your csv files separating data from column headers/text
    % if you get an 'unable to open file error, make sure that the
    % folder/files are on the matlab path
    if length(raw_data(:,:,File).data) < (aq_freq*tr_length-1) %makes sure file is correct length
         DLC_data(:,:,File) = zeros((aq_freq*tr_length),16); % if file is incorrect length pads matrix with zeros
    else
    DLC_data(:,:,File) = (raw_data(:,:,File).data(:,[ 2 3 5 6 8 9 11 12 14 15 17 18 20 21 23 24])); 
    % grabs just the data you want from each file 
    %(x,y position of eye) and puts it into a 3D matrix where the first two dimensions are the data and the third dimension is which file it came from
    end
end
if length(DLC_data(:,1,:)) > 6600; %figures out if videos are 20Hz or 40Hz automatically by frames cutoff
    fs = 40
else
    fs = 20; % video sampling rate in Hz
end
    time = ((1:length(DLC_data))*(1/fs)+1/fs)-(1/fs); % creates a time vector according to the sampling rate and length of video starting at 0.03s
    time = time'; % transposes time vector from row to column vector
%%  this part of the code fits a circle to the point centers for each frame
    circData = zeros(length(DLC_data),3,size(DLC_data,3)); % initializes the circData variable
    points = zeros(8,2); % initializes the points variable that will contain the centroid of each track point
    pupil_max=[];
    for z = 1:size(DLC_data,3);
        TF(:,:,z)=isoutlier(DLC_data(:,:,z),'gesd');
            for i = 1:length(TF);
                for j = 1:size(TF,2);
                    if TF(i,j,z)==1;
                        DLC_data(i,j,z) = NaN;
                    end
                end
            end
    end
    for x = 1:size(DLC_data,3); % moves in the 3rd dimension (trial number in this case)
        for i = 1:size(DLC_data,1) % moves along each video frame
            for j = 1:8 %moves along the centroid of each tracking point
                points(j,1,x) = DLC_data(i,(j*2)-1,x); % x value for the centroid of track points
                points(j,2,x) = DLC_data(i,(j*2),x); % y value for the centroid of each tracked point
            end
            circData(i,:,x) = CircleFitByPratt(points(:,:,x)); % fit circle to tracked points at each video frame...
            % and for each trial
        end
        pupil_diameter(:,:,x)=horzcat([time (circData(:,3,x)*2)]);
%         for ix = 2:length(pupil_diameter);
%             if pupil_diameter(:,2,x)==NaN;
%                 a = find(pupil_diameter(ix,2,x),1,'last')
%                 pupil_diameter(ix,2,x) = pupil_diameter
%             else
%                 pupil_diameter(ix,2,x)=pupil_diameter(ix,2,x);
%             end
%         end
        m = max(pupil_diameter(:,2,x));
        pupil_max = [pupil_max m];
    end 
    
%% save the data
dir_out = uigetdir(); % asks you where you want to save the new files
save([dir_out,'\',sprintf('%06s_%11s_pupil_track_filt',date,mouseID)],'DLC_data','dir_in','circData','pupil_diameter','pupil_max'); %saves DLC data in matlab format
% saves in chosen folder with date and mouse ID

%% let's plot the data!
x = input('Which trials do you want to plot?');
for i = 1:length(x);
    plot(pupil_diameter(:,1,x(i)),(pupil_diameter(:,2,x(i))/pupil_max(x(i)))*100,'LineWidth',1);
    hold on;
end

% if strcmp(dir_in(31:34),'Halo');
%     barcolor=[0.8 0.3 0];
% else
%     barcolor=[0 0 1];
% end
% stim_start = input('enter stim start time');
% stim_end = input('enter stim end time');
% a = [stim_start stim_end stim_end stim_start];
% b =[min(ylim) min(ylim) max(ylim) max(ylim)]; 
% fill(a,b,barcolor,'FaceAlpha',0.5);
set(gca,'LineWidth',2,'TickDir','out','TickLength',[0.025 0.025],'FontSize',14);
box off;
xlabel('time (s)');
ylabel('pupil diameter (% max)');
xlim([0 330]);
