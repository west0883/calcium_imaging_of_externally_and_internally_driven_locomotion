

function [parameters] = FindBetaChanges(parameters)

    MessageToUser('Beta change of ', parameters);

    % Inputs
    dataset = parameters.dataset; % Output from DatasetPrep

    if ~parameters.onPermutations
        if strcmp(parameters.comparison_type, 'categorical')
            betas = parameters.betas(2, :); % row 1 is intercepts; 
        else
            betas = parameters.betas(2:end, :); % row 1 is intercepts; 
        end 
    else 
        if strcmp(parameters.comparison_type, 'categorical')
            betas =  parameters.betas;
        else 
            betas = parameters.betas(2:end, :, :); % row 1 is intercepts; 
        end 
    end 


    % Inverted
    ysig = dataset.zscoring.explanatoryVariables.sigma;

    if strcmp(parameters.comparison_type, 'categorical')
        xsig = dataset.zscoring.responseVariables.sigma(1); % take only first because categorical
    else 
        xsig = dataset.zscoring.responseVariables.sigma;
    end 

    % Make dimensions match (replicate xsig so there's a set for each response varaible.
    xsig = repmat(xsig, size(ysig,2),1); 
  

   % Calculate sigmas
   sigmas = ysig./ xsig';

   % multiply betas by sigmas
   beta_change = betas .* sigmas;

   parameters.beta_change = beta_change;

end 