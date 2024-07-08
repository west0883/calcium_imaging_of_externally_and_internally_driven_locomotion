% AverageByNode.m
% Sarah West
% 7/25/22

% Takes correlation-type data and finds an average of that data by node.
% Uses RunAnalysis.
% Inputs: 
% parameters.data -- correlation values (correlation matrix or lower triangle of a
% correlation matrix only). Matrices must be put in with number_of_sources
% x number_of_sources dimensions as first two dimensions.
% parameters.isVector -- flag saying if the inputted data is in the
% lower-triangle vector form or not. If it is, then the vector needs to be
% reshaped into a matrix to identify which correlations belong to which
% node pairs.
% parameters.number_of_sources -- number of unique nodes.

function [parameters] = AverageByNode_Inverted(parameters)

    % Send message to user.
    MessageToUser('Averaging ', parameters);

    number_of_sources = parameters.number_of_sources;

    % Find if this dataset was from PLSR (and needs to be separated by
    % variable type)
    if isfield(parameters, 'fromPLSR') && parameters.fromPLSR
        PLSR_flag = true;

    else
        PLSR_flag = false;

    end

    % Calculate indices of the lower triangle if hasn't been done already.
    if ~isfield(parameters, 'indices')
        parameters.indices = find(tril(ones(parameters.number_of_sources), -1));
    end
  
    % Pull out data from parameters structure for easier use.
    data = parameters.data;

    % If using average sigmas, pull those out, too
    if parameters.multiply_by_average_sigma
        sigmas = parameters.average_sigmas;
    end

    % If correlations are in vector format, reshape into a full correlation
    % matrix.
    if parameters.isVector

        % If the different correlation pairs aren't in the first dimension,
        % permute so they are.
        if parameters.corrsDim ~= 1
            
            data = permute(data, [parameters.corrsDim  setdiff(1:ndims(data), parameters.corrsDim, 'stable')]);
           
        end 

        if parameters.multiply_by_average_sigma
            if parameters.sigmasDim ~= 1
                
                sigmas= permute(sigmas, [parameters.corrsDim  setdiff(1:ndims(data), parameters.corrsDim, 'stable')]);
               
            end 
        end

        % If user says to only use the changes found to be significant, replace
        % non-significant values with NaNs
        if isfield(parameters, 'significantOnly') && parameters.significantOnly
    
            if parameters.isVector
                data(~parameters.significance,:,:) = NaN;
            else
    
                error('no code yet for significantOnly with non-vector formats');
    
            end
    
        end


        % If this is a "comparison" from PLSR, need to check if there was
        % more than one variable included (like for "continuous"
        % comparisons).
        if PLSR_flag
       

            % *** REMOVE THIS STEP FOR INVERTED *** 
            % If Betas, remove intercepts from data
%             if strcmp(output_type, 'BETA') 
%     
%                C = repmat({':'}, 1, ndims(data));
%                C{1} = [1:numel(parameters.indices) + 1:size(data, 1)];
%                data(C{:}) = [];
%             
%             end 
            % Find number of response variables (will be a multiple of
            % number of unique correlation pairs).

            unique_corr_pairs = number_of_sources * (number_of_sources - 1)/2;
            response_variable_number = size(data, 1) / unique_corr_pairs;

            if rem(size(data,1), unique_corr_pairs) ~= 0
                error('Wrong number of correlations.')
            end
            
            data = reshape(data, unique_corr_pairs, response_variable_number, size(data,2), size(data,3));
            if parameters.multiply_by_average_sigma
                sigmas =  reshape(sigmas, unique_corr_pairs, response_variable_number, size(sigmas,2), size(sigmas,3));
            end

        end

        % Multiply data by sigmas
        if parameters.multiply_by_average_sigma
                data = data .* sigmas;
        end


        % Make a holder matrix, insert the data.
        % For now just assuming the original vectors had no more than 3
        % dimensions.
        holder = NaN(parameters.number_of_sources, parameters.number_of_sources, size(data, 2), size(data, 3), size(data,4));

        % Get some variables before parfor loop
        indices = parameters.indices;
        size3 = size(holder,3);
        size4 = size(holder, 4);
        size5 = size(holder, 5);

        indices_upper = find(triu(ones(parameters.number_of_sources), 1));

        % Insert.
        for k = 1:size5
            for i = 1:size3
                for j = 1:size4

                    subholder = NaN(number_of_sources, number_of_sources);
                
                    subholder(indices) = data(:, i, j,k);
                   
    
                    % Duplicate betas across diagonal.
                    betas_flipped = subholder';
                    elements_upper = betas_flipped(indices_upper);
                    subholder(indices_upper) = elements_upper;
    
                    holder(:, :, i, j, k) = subholder;
                end
            end
        end

        % Rename holder to match non-vector condition.
        data = holder;

    else
       
        % If not a vector, make sure there aren't any values in the
        % diagonals.
        for i = 1:number_of_sources
            data(i, i, :, :) = NaN;
        end
    end

    % Check if there are any infinities in the matrix (can happen if values
    % of 1 are sent through the Fisher transform). Turn to NaNs if there are. 
    index = data == Inf;
    if any(index, 'all')
        data(index) = NaN;
    end

    % Take averages.Also calculate sums (is useful when you don't want non-significant
    % values.)
    node_averages = mean(data, 1, 'omitnan');
    node_stds = std(data, [], 1, 'omitnan');
    node_sums = sum(data, 1, 'omitnan');

    % If from PLSR, return different variables to same dimension (will
    % probably work better with later steps that were originally designed
    % this way). 
    if PLSR_flag

        node_averages = reshape(node_averages, [1, number_of_sources * response_variable_number, size(node_averages,4), size(node_averages,5)]);
        node_stds = reshape(node_stds, [1, number_of_sources * response_variable_number, size(node_stds,4), size(node_stds,5)]);
        node_sums = reshape(node_sums, [1, number_of_sources * response_variable_number, size(node_sums,4), size(node_sums,5)]);
    end

    % Now permute out the leading 1.
    node_averages = squeeze(node_averages);
    node_stds = squeeze(node_stds);
    node_sums = squeeze(node_sums);

    % Put into output structure.
    parameters.node_averages = node_averages;
    parameters.node_stds = node_stds;
    parameters.node_sums = node_sums;

end