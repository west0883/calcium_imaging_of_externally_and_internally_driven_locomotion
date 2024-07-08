% IpsaContraAverage
% Sarah West
% 7/22/22
% Takes a correlation matrix and averages the values of nodes with their 
% corresponding node on the opposite hemisphere, but keeping ipsalateral & 
% contralateral relationships intact. Reduces the number of unique
% correlations by half.

function [parameters] = IpsaContraAverage(parameters)
    
    MessageToUser('Ipsa-contra on ', parameters);

    % If empty, skip this value.
    if isempty(parameters.data)
        parameters.data_ipsacontra = [];
        return
    end

    holder = NaN(size(parameters.data));
    lefts  = 1:2:parameters.number_of_sources;
    rights = 2:2:parameters.number_of_sources;
    
    % ipsa
    ipsas_left = parameters.data(lefts,1:2:parameters.number_of_sources, :,:);
    ipsas_right = parameters.data(rights, 2:2:parameters.number_of_sources, :, :);
    
    % contra/even
    contras_left = parameters.data(lefts, 2:2:parameters.number_of_sources, :, :);
    contras_right = parameters.data(rights, 1:2:parameters.number_of_sources, :,:);
    
    % Take means
    contras = mean(cat(5, contras_left, contras_right), 5);
    ipsas = mean(cat(5,ipsas_left, ipsas_right), 5);
    
    % put back in
    % ispas, left then right
    holder(lefts, 1:2:parameters.number_of_sources, :,:) = ipsas;
    holder(rights, 2:2:parameters.number_of_sources,:,:) = ipsas;
    
    % contras, left then right
    holder(lefts, 2:2:parameters.number_of_sources, :, :) = contras;
    holder(rights, 1:2:parameters.number_of_sources, :,:) = contras;
    
    % Put in the correlations between the homologous regions
    % holder(rights, lefts) = parameters.data(rights, lefts);
    for i = 1:numel(rights)
       holder(rights(i), lefts(i),:,:) = parameters.data(rights(i), lefts(i),:,:);
    end

    % Put into output matrix, make single precision for saving.
    parameters.data_ipsacontra = single(holder);

end
