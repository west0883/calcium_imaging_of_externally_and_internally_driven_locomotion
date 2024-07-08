% DivideData.m
% Sarah West
% 5/10/22

% Inputs: 
% parameters.division_points -- an array of index points along divideDim where
% the data should be divided.
% parameters.divideDim -- a scalar. The dimension you're dividing across. 
function [parameters] = DivideData(parameters)

    % Display progress message to user.
    MessageToUser('Dividing ', parameters);

    % If divisions were from concatenate data, (will be from
    % concatenated_origin variables)
    if isfield(parameters, 'fromConcatenateData') && parameters.fromConcatenateData
       
        % Write this in multiple lines because that's how cell arrays work.
        holder = [parameters.division_points{:}];
        division_points = cumsum([holder{1,:}]);
    
    % If not from ConcatenateData (default)
    else 
        division_points = parameters.division_points;
    end 

    % If division points doesn't have a 0 as the first point, add a 0 to
    % beginning. 
    if division_points(1) ~= 0

        % Add a zero to beginning, even if row or column vector. 
        division_points = [0 division_points(1:end)];

    end
 
    % Make a holder for divided data (don't count the 0)
    data_divided = cell(numel(division_points) - 1, 1); 

    % Make flexible dimensionts. 
    C = repmat({':'},1, ndims(parameters.data));

    % Divide the data, putting each division into its own entry of a cell array.
    % Skip first entry because it's 0.
    for dividei = 2:numel(division_points)

        C(parameters.divideDim) = {division_points(dividei - 1) + 1 : division_points(dividei)};
        data_divided{dividei - 1} = parameters.data(C{:});
    end  

    % Put into output 
    parameters.data_divided = data_divided;
end 