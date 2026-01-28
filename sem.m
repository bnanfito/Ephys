function [out] = sem(x)
    % input x should be either a vector representing observations of a
    % single variable or a 2D matrix of observations for multiple variables
    
    % if x is a matrix, each row is treated as an observation and each
    % column as a variable
    if size(x,1)>1 && size(x,2)>1
        sig = std(x,[],1,'omitnan');
        n = size(x,1);
    else
        sig = std(x,'omitnan');
        n = length(x);
    end
    out = sig/sqrt(n);

end