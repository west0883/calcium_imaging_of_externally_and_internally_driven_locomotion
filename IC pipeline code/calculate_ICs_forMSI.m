% calculate_ICs.m
% Sarah West
% 9/1/21

% Calculates ICs from SVD compressed data. Assumes one compressed dataset
% per mouse.

function []= calculate_ICs_forMSI(mouse_number) 
    
    % Parameters.
    mouse = num2str(mouse_number);
    num_sources = 100;
    spatial_component = 'V';
     
    folder=pwd;
    %addpath(genpath(folder));

    dir_in=[folder '/']; % directory on the MSI network. 
    filename_input = [dir_in 'm' mouse '_SVD_compressed.mat'];
    dir_out=[folder '/' ]; % directory on the MSI network.
    filename_output = [dir_out 'm' mouse '_sources' num2str(num_sources) '.mat']; 

    % Tell us the mouse number
    disp(['mouse #' mouse]); 
        
    % Load in the compressed data-- only the spatial component and the S. 
    load(filename_input, 'S', spatial_component);
    
    % Depending on which component is the spatial component, calculate 
    % the sources accordingly. 
    switch spatial_component
        case 'V'
            B=jader_lsp([S*V'],num_sources);
            sources=B*[S*V'];
        case 'U'
            B=jader_lsp([U*S],num_sources);
            sources=B*[U*S];
    end
    
    % Create output file path & filename
    mkdir(dir_out); 

    % Save sources and B.
    save(filename_output, 'sources', 'B', '-v7.3');  
end