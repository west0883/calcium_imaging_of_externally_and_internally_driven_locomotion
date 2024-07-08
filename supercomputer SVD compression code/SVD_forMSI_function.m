% SVD_forMSI_function.m
% Sarah West 
% 9/6/21

% Runs SVD as a function with mouse number as an input

% Runs SVD on all the data from a given animal. Must be done with MSI or everything will die.
% For each mouse, loads the filtered and blood vessel regressed data, concatenates it all
% together, runs SVD, then saves the outputs.

function []=SVD_forMSI_function(mouse_number)
    
    % Convert mouse number to a string 
    mouse=num2str(mouse_number); 
    
    n_compressions = 500;
    
    folder=pwd;
    %addpath(genpath(folder));
    dir_in=[folder '/fully preprocessed stacks/' mouse '/']; % directory on the MSI network. 
    dir_out=[folder '/' ]; % directory on the MSI network.
    
    % Load the list of days included for each mouse.
    load([folder '/mice_all_random.mat']); 

    % Determine index of mouse within mice_all.
    mousei = find(strcmp({mice_all.name}, mouse)==1);
    
    disp(['mouse ' mouse]);

    % Make output filename
    filename_output=[dir_out 'm' mouse '_SVD_compressed.mat']; 

    % Start a running count of how many stacks each mouse has; for data space pre-alotment
    total_stacks=0;
    stacks = {};
    % For each  day; count how many stacks are in each day and add them all up so you can make an accurately sized matrix for data pre-alotment
    for dayi=1:size(mice_all(mousei).days,2)  
        
        % Get the day name.
        day=mice_all(mousei).days(dayi).name; 
        
        all_stacks = [mice_all(mousei).days(dayi).stacks mice_all(mousei).days(dayi).spontaneous];
        all_stacks = all_stacks(~isnan(all_stacks));
        
        % Get list of stacks for that day from mice_all
        for stacki = 1: numel(all_stacks)
            % List the stacks in a given day
            filename = sprintf('data%02d.mat', all_stacks(stacki));
            stacks= [stacks; [dir_in day '/' filename]]; 
        end

        % Add the number of stacks to the running count
        total_stacks=total_stacks+size(all_stacks,2); 
    end 

    disp(['total stacks =' num2str(total_stacks)]); 
    
    % Get infor from first stack to get nummber of pixels
    day=mice_all(mousei).days(1).name; 
     
    matObj = matfile(stacks{1});
    pixels = size(matObj,'data', 1);
    frames=size(matObj, 'data' ,2);
   
    % Initialize data matrix.
    all_data=NaN(total_stacks*frames, pixels); 

    % load and concatenate data
    disp('concatenating');
    count=0; 
    
    % for each stack
    for stacki=1:size(stacks,1) 
        load(stacks{stacki});  
        all_data((count*frames)+1:(count+1)*frames, :)=data';
        count=count+1;
    end 

    % run SVD   
    disp(['running SVD, ' num2str(n_compressions) ' compressions']); 
    [U,S,V]=compute_svd(all_data, 'randomized', n_compressions);
    disp('saving'); 
    save(filename_output, 'U', 'S', 'V', '-v7.3'); 

end
