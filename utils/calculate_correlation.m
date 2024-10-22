function CC = calculate_correlation(Target, Input, mask_type, flag, nonzero, corr_type)
    % Input:
    % Target: The target matrix (square matrix)
    % FC: The functional connectivity matrix (square matrix with the same dimensions as Target)
    % mask_type: The mask type (1 for original mask, 2 for full matrix without diagonal)
    % flag: If set to 1, a scatter plot of the input and target variables with the best fit line will be created
    % Output:
    % CC: The adjusted R-squared value representing the correlation between the target and functional connectivity matrices

    % Check if the input matrices have the same dimensions
    if ~isequal(size(Target), size(Input))
        error('The target and functional connectivity matrices must have the same dimensions.');
    end
    

    % Calculate the mask
    dim = size(Target, 2);
    if mask_type == 1
        mask = triu(ones(dim, dim), 0);
    elseif mask_type == 2
        mask = ~eye(dim);
    else
        mask = ones(length(Target),1);
    end

    % Calculate the correlation between the target and functional connectivity matrices
    if nonzero == 1
        indices = (Target ~=0) & (Input ~= 0) & mask;
    elseif nonzero == 2
        indices = (Target ~=0) & mask;
    else
        indices = mask;
    end
    
    if size(Target, 1) == size(Target, 2)
        x1 = Input(indices);
        y1 = Target(indices);
    else
        x1 = Input;
        y1 = Target;
    end
    mdl = fitlm(x1, y1);
    if strcmp(corr_type,'Rsquared')
        CC = mdl.Rsquared.adjusted;
    elseif strcmp(corr_type,'Spearman')
        CC = corr(x1,y1,'Type','Spearman');
    elseif strcmp(corr_type,'Pearson')
        CC = corr(x1,y1);
    end

    % Create the scatter plot if flag is set to 1
    if flag == 1
        figure;
        scatter(x1, y1,'SizeData',10);
        hold on;
        plot(x1, mdl.Fitted, 'r--');
        hold off;
        xlabel('Input');
        ylabel('Target');
        title('Scatter plot with best fit line');
    end

end
