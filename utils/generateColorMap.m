function [cMap] = generateColorMap(data, resolution)
    % generateColorMap - Generate a color map based on data values.
    % This function calculates a color map that transitions from blue to red
    % based on the minimum and maximum values of the input data.
    %
    % Syntax: [cMap] = generateColorMap(data, resolution)
    %
    % Inputs:
    %    data - A matrix containing the data values.
    %    resolution - The number of colors to be used in the color map.
    %
    % Outputs:
    %    cMap - A matrix representing the generated color map.
    
    % Calculate the minimum and maximum values from the data
    cmin = min(data(:));
    cmax = max(data(:));

    % Initialize the number of blue and red colors
    nr_blues = 0;
    nr_reds = 0;

    % Calculate the sum of the absolute values of cmin and cmax for normalization
    sum_abs_values = abs(cmin) + abs(cmax);

    % Determine the number of blue and red colors based on the data range
    if cmin < 6
        if cmax > 0 % cmax is positive
            % Normalize and calculate the number of blues and reds
            nr_blues = round(resolution * (abs(cmin) / sum_abs_values));
            nr_reds  = round(resolution * (abs(cmax) / sum_abs_values));
        else % cmax is also negative
            nr_blues = round(resolution);
            nr_reds  = 0;
        end
    else
        % cmin is non-negative (and thus, cmax is also non-negative)
        nr_blues = 0;
        nr_reds  = round(resolution * (abs(cmax) / max(cmin, cmax)));
    end

    % Generate the color map based on the calculated number of blues and reds
    if nr_blues == 0
        cMap = red2black(nr_reds);
        % cMap = flipud(red2white(nr_reds)); % Uncomment to use red to white transition
    elseif nr_reds == 0
        cMap = flipud(white2blue(nr_blues));
    else
        cMap = [cyan2black(nr_blues); black2yellow(nr_reds)];
        % cMap = [flipud(white2blue(nr_blues)); flipud(red2white(nr_reds))]; % Uncomment to use white to blue and red to white transition
    end
end
