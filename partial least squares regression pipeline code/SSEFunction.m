% SSEFunction.m
% Sarah West
% 6/30/22

% Function that calculates sum of squares on a PLSR model. Taken from
% source code of Matlab built-in function plsregress.

function sumsqerr = SSEFunction(Xtrain,Ytrain,Xtest,Ytest,ncomp)

    XmeanTrain = mean(Xtrain);
    YmeanTrain = mean(Ytrain);
    X0train = Xtrain - XmeanTrain;
    Y0train = Ytrain - YmeanTrain;
    
    % Get and center the test data
    X0test = Xtest - XmeanTrain;
    Y0test = Ytest - YmeanTrain;
    
    % Fit the full model, models with 1:(ncomp-1) components are nested within
    [Xloadings,Yloadings,~,~,Weights] = SIMPLSFunction(X0train,Y0train,ncomp);
    XscoresTest = X0test * Weights;
    
    % Return error for as many components as the asked for.
    outClass = superiorfloat(Xtrain,Ytrain);
    sumsqerr = zeros(2,ncomp+1,outClass); % this will get reshaped to a row by CROSSVAL
    
    % Sum of squared errors for the null model
    sumsqerr(1,1) = sum(sum(abs(X0test).^2, 2));
    sumsqerr(2,1) = sum(sum(abs(Y0test).^2, 2));
    
    % Compute sum of squared errors for models with 1:ncomp components
    for i = 1:ncomp
        X0reconstructed = XscoresTest(:,1:i) * Xloadings(:,1:i)';
        sumsqerr(1,i+1) = sum(sum(abs(X0test - X0reconstructed).^2, 2));
    
        Y0reconstructed = XscoresTest(:,1:i) * Yloadings(:,1:i)';
        sumsqerr(2,i+1) = sum(sum(abs(Y0test - Y0reconstructed).^2, 2));
    
    end
end