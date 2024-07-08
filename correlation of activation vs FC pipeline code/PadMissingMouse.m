% PadMissingMouse.m
% Sarah West
% 1/20/23

% A quick and dirty function to concatenate data across conditions for the 
% correlation fluor vs corr pipeline

function [parameters] = PadMissingMouse(parameters)

    % Inputs :
    % parameters.correlations -- 7 x 16 (mice x unique nodes)
    % parameters.fluorescence -- 7 x 16 (mice x unique nodes)
    % parameters.mouseDim -- scalar; dimension of correlations and
    % fluorescence the different mice are in.
    % parameters.placement -- scalar; where in the order hte missing mouse
    % should go

    % take these out becuase you might be altering them
    correlations = parameters.correlations;
    fluorescence = parameters.fluorescence;
    mouseDim = parameters.mouseDim;


    % if mouse 1100 isn't there, 
    if size(correlations, mouseDim) ~= parameters.number_of_mice

        sizes = size(correlations);

        % Create matrices for variying dimension indexing. 
        Cin= repmat({':'}, 1, numel(sizes));
        Cout = repmat({':'}, 1, numel(sizes));

        sizes = size(correlations);
        sizes_new = sizes;
        sizes_new(mouseDim) = sizes(mouseDim) + 1;

        correlations_padded  = NaN(sizes_new);
        fluorescence_padded = NaN(sizes_new);

        % Put in data
        Cout{mouseDim} = [1:parameters.placement - 1, parameters.placement + 1:sizes_new(mouseDim)]; 
      
        correlations_padded(Cout{:}) = correlations; 
        fluorescence_padded(Cout{:}) = fluorescence;
        
    else 
        correlations_padded = correlations;
        fluorescence_padded = fluorescence;
    end 

    parameters.corrs_padded = correlations_padded;
    parameters.fluors_padded = fluorescence_padded;

end 