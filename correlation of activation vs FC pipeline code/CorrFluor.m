% CorrFluor.m
% Sarah West
% 1/20/23

% High-level function to calculate correlation between changes in 
% fluorescence and correlation per condition. 

function [parameters] = CorrFluor(parameters)

    % parameters.correlations -- 7 x 16 x 48 ; mouse, node, comparison
    % parameters.fluorescence -- 7 x 16 x 48 ;   mouse, node, comparison\
    % parameters.mouseDim 
  
    mouseDim =parameters.mouseDim;
    nodeDim = parameters.nodeDim;
   
    correlations = parameters.correlations;
    fluorescence = parameters.fluorescence;

    sizes = size(correlations);

    corrs_per_mouse = NaN(sizes(1:end-1));
    corrs_per_mouse_pvalues = NaN(sizes(1:end-1));

    % If not continuous 
    if ~ strcmp(parameters.comparison_type, 'continuous')
   
        % for each mouse
        for mousei = 1:sizes(parameters.mouseDim)
            
            % for each node
            for nodei = 1:sizes(parameters.nodeDim)

                Cin= repmat({':'}, 1, numel(sizes));
                Cout = repmat({':'}, 1, numel(sizes) - 1);

                Cin{mouseDim} = mousei;
                Cin{nodeDim} = nodei; 
                Cout{mouseDim} = mousei;
                Cout{nodeDim} = nodei; 
    
                fluor = squeeze(fluorescence(Cin{:}));
                corr = squeeze(correlations(Cin{:}));
              
                % remove any nans
                fluor(isnan(fluor)) = [];
                 
                corr(isnan(corr)) = [];

                % Run the correlation
    
                [r, p] = corrcoef(fluor, corr);
    
                corrs_per_mouse(Cout{:}) = r(1, 2); 
                corrs_per_mouse_pvalues(Cout{:}) = p(1, 2); 
    
            end 
    
        end

    % If continuous,
    else
        variableDim = parameters.variableDim;
        
        % for each mouse
        for mousei = 1:sizes(parameters.mouseDim)
            
            % for each varible 
            for variablei = 1:sizes(variableDim)

                % for each node
                for nodei = 1:sizes(parameters.nodeDim)
    
                    Cin= repmat({':'}, 1, numel(sizes));
                    Cout = repmat({':'}, 1, numel(sizes) - 1);
    
                    Cin{mouseDim} = mousei;
                    Cin{variableDim} = variablei; 
                    Cin{nodeDim} = nodei; 
                    Cout{mouseDim} = mousei;
                    Cout{variableDim} = variablei;
                    Cout{nodeDim} = nodei; 
                   
                    fluor = squeeze(fluorescence(Cin{:}));
                    corr = squeeze(correlations(Cin{:}));
                  
                    % remove any nans
                    fluor(isnan(fluor)) = [];
                     
                    corr(isnan(corr)) = [];
    
                    % Run the correlation
        
                    [r, p] = corrcoef(fluor, corr);
        
                    corrs_per_mouse(Cout{:}) = r(1, 2); 
                    corrs_per_mouse_pvalues(Cout{:}) = p(1, 2); 
                end
            end 
        end
    end 
    % Put into output structure
    parameters.corrs_per_mouse = corrs_per_mouse;
    parameters.corrs_per_mouse_pvalues = corrs_per_mouse_pvalues;

    % Average across mice
    parameters.corrs_across_mice = squeeze(mean(corrs_per_mouse, 1, 'omitnan'));
   

    
end 