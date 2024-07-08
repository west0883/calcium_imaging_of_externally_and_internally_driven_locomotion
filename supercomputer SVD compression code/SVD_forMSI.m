% SVD_forMSI.m
% Sarah West 
% 8/6/20

% Runs SVD on all the data from a given animal. Must be done with MSI or everything will die.
% For each mouse, loads the filtered and blood vessel regressed data, concatenates it all
% together, runs SVD, then saves the outputs.
clear all;
n_compressions=200;
pixels=256; %size of images in your stacks in pixels, in one dimension (assuming the stack is square)
frames=6000; % number of frames in each stack 

folder=pwd;
addpath(genpath(folder));
dir_in=[folder '/' ]; % directory on the MSI network. 
dir_out=[folder '/' ]; % directory on the MSI network.

% Load list of days and mice. 
load([dir_in 'days_all.mat']);

% For each mouse
for mousei=2 %:3 %4:size(mice,1)   
    mouse=mice{mousei};
    disp(['mouse #' mouse]);
    eval(['days_list=days_all.m' mouse ';']); 
    total_stacks=0; % start a running count of how many stacks each mouse has; for data space pre-alotment
    for dayi=1:size(days_list,1) % for each  day; count how many stacks are in each day and add them all up so you can make an accurately sized matrix for data pre-alotment
        day=days_list(dayi,:);
        stacks=dir([dir_in 'preprocessed round 4/' day '/data*.mat']);  % list the  stacks in a given day
        total_stacks=total_stacks+size(stacks,1); % add the number of stacks to the running count
    end 
    disp(['total stacks =' num2str(total_stacks)]); 
    all_data=NaN(total_stacks*frames, pixels*pixels);  % initialize data matrix
    
    load and concatenate data
    disp('concatenating');
    count=0;
    try
        for dayi=1:size(days_list,1) % for each  day; 
            day=days_list(dayi,:);
            stacks=dir([dir_in 'preprocessed round 4/' day '/data*.mat']);  % list the  stacks in a given day
            for stacki=1:size(stacks,1) % for each stack
                load([dir_in 'preprocessed round 4/' day '/' stacks(stacki).name]);  
                data_reshaped=reshape(data,256*256, 6000); 
                all_data((count*6000)+1:(count+1)*6000, :)=data_reshaped';
                count=count+1;
            end 
        end
    % run SVD   
    disp('running SVD'); 
    [U,S,V]=compute_svd(all_data, 'randomized', n_compressions);
    save([dir_out 'm' mouse '_SVD_compressed_round4processing.mat'], 'U', 'S', 'V', '-v7.3'); 
    catch
        disp('found a corrupt file');
    end
    clearvars -except n_compressions dir_in dir_out mice mousei days_all
end
