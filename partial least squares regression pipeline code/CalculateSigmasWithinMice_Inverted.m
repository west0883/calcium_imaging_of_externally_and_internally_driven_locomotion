% CalculateSigmasWithinMice_Inverted.m
% Sarah West
% 10/24/23

% need for calculating corr of fluorescence with corrs by nodes.
% Made to be called by RunAnalysis.m
function [parameters] = CalculateSigmasWithinMice_Inverted(parameters)


    MessageToUser('Calculating sigmas ', parameters);

        % INVERT HERE 
       ysig = parameters.dataset.zscoring.explanatoryVariables.sigma;


       % Get x sigma
       % If using just the x zscore, make xsigma = 1.
       % INVERT HERE

       % Keep only the first response variable if comparison type is
       % categorical.
       xsig_single = parameters.dataset.zscoring.responseVariables.sigma;

       % get comparison type 
       comparison_type = parameters.values{strcmp(parameters.keywords, 'comparison_type')};
       
       if strcmp(comparison_type, 'categorical')

           xsig_single = xsig_single(1);

       end

       % Make dimensions match (replicate xsig so there's a set for each response varaible.
       xsig = repmat(xsig_single, size(ysig,2),1); 

       % Calculate sigmas
       sigmas = reshape((ysig ./xsig')', 1, []);

       % Put into output structure
       parameters.sigmas = sigmas;
end 