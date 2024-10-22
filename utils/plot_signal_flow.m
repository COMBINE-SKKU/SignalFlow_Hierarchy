function plot_signal_flow(nModules, signal_flow_pos, signal_flow_neg)
% PLOT_SIGNAL_FLOW Plots the positive and negative signal flow for a given number of modules.
%
%   plot_signal_flow(nModules, signal_flow_pos, signal_flow_neg)
%
%   Inputs:
%       nModules        - The number of modules to plot.
%       signal_flow_pos - A vector containing the positive signal flow for each module.
%       signal_flow_neg - A vector containing the negative signal flow for each module.
%
%   This function generates a plot showing the positive and negative signal flows
%   between modules, with edge bundling applied to enhance visualization.

    % Generate X and Y coordinates for plotting
    Xcorrds = linspace(-150, 150, nModules);  % X coordinates from -150 to 150 with nModules points
    Ycorrds = [ones(1, nModules) * 50, ones(1, nModules) * -50];  % Y coordinates alternating between 50 and -50

    % Replicate X coordinates for two sets (for modules and their targets)
    x = repmat(Xcorrds, 1, 2);  % Xcorrds repeated twice for creating a full X vector
    y = Ycorrds;  % Y coordinates for plotting

    % Create IDs for modules and their corresponding targets
    id = arrayfun(@(n) sprintf('%d', n), 1:nModules, 'UniformOutput', false);  % IDs without 't'
    id_t = strcat(id, 't');  % Append 't' to IDs to create target IDs
    id_all = [id, id_t];  % Combine both lists

    % Generate source (src) and target (tar) strings for the graph edges
    src = repelem(id, nModules);  % Repeat each ID nModules times
    tar = repmat(id_t, 1, nModules);  % Repeat the entire target list nModules times

    % ----------------------------- Plot Negative Edges ----------------------------- %
    % Calculate the absolute weights of negative signal flow and create a directed graph
    w_neg = abs(signal_flow_neg(:))';  % Get the absolute values of signal flow (negative edges)
    G_neg = digraph(src, tar, w_neg);  % Create a directed graph with source, target, and weights

    % Map node positions to graph nodes based on the id order
    [~, loc] = ismember(G_neg.Nodes.Name, id_all);  % Find the locations of node names in the ID array
    G_neg.Nodes = [G_neg.Nodes table(x(loc)', y(loc)', 'VariableNames', {'x', 'y'})];  % Add x and y coordinates to nodes

    % Plot the negative signal flow graph
    figure;
    h_neg = plot(G_neg, 'XData', G_neg.Nodes.x, 'YData', G_neg.Nodes.y, ...
        'MarkerSize', 10, 'NodeLabel', []);  % Plot the graph without node labels
    axis equal;  % Ensure the axes are scaled equally

    % Enhance edge visualization for negative signal flow using edge bundling
    he2_neg = plotdeb(G_neg, 'initial', true);  % Plot the initial debundled edges
    set(h_neg, 'EdgeColor', 'none');  % Remove edge colors from the original plot
    uistack(h_neg, 'top');  % Bring the node plot to the top

    % Debundle the edges to reduce overlap and enhance visibility
    G_neg = debundle(G_neg);  % Debundle the graph to reduce edge overlaps
    delete(he2_neg);  % Delete the initial debundled edges
    he3_neg = plotdeb(G_neg, 'w', 2);  % Plot the debundled edges with a specified width
    uistack(h_neg, 'top');  % Bring the node plot to the top

    % Customize the appearance of the negative edges
    set(he3_neg, 'FaceColor', "#2121cf");  % Set the edge color to a specific blue shade
    box off;  % Remove the box around the plot
    set(gca, 'Visible', 'off');  % Hide the axes

    % ----------------------------- Plot Positive Edges ----------------------------- %
    % Calculate the weights of positive signal flow and create a directed graph
    w_pos = signal_flow_pos(:)';  % Get the signal flow (positive edges)
    G_pos = digraph(src, tar, w_pos);  % Create a directed graph with source, target, and weights

    % Map node positions to graph nodes based on the id order
    [~, loc] = ismember(G_pos.Nodes.Name, id_all);  % Find the locations of node names in the ID array
    G_pos.Nodes = [G_pos.Nodes table(x(loc)', y(loc)', 'VariableNames', {'x', 'y'})];  % Add x and y coordinates to nodes

    switch nModules
        case 22
            % Load custom colors for the nodes 
            clrs = MMP22moduleColors;  
            colors = [clrs(2,:); clrs(2:end,:); clrs(3:end,:)];  
        case 27
             % Load custom colors for the nodes 
            clrs = importdata('/combinelab/03_user/younghyun/01_project/01_HierarchyMapping/data/module27(ordered)clrs.mat');
            colors = [clrs(1,:);clrs(1:end,:);clrs(2:end,:)];
    end
            
    % Plot the positive signal flow graph with custom node colors
    h_pos = plot(G_pos, 'XData', G_pos.Nodes.x, 'YData', G_pos.Nodes.y, ...
        'MarkerSize', 10, 'NodeLabel', [], 'NodeColor', colors);  

    % Enhance edge visualization for positive signal flow using edge bundling
    he1_pos = plotdeb(G_pos, 'initial', true);  % Plot the initial debundled edges
    set(h_pos, 'EdgeColor', 'none');  % Remove edge colors from the original plot
    uistack(h_pos, 'top');  % Bring the node plot to the top

    % Debundle the edges to reduce overlap and enhance visibility
    G_pos = debundle(G_pos);  % Debundle the graph to reduce edge overlaps
    delete(he1_pos);  % Delete the initial debundled edges
    he2_pos = plotdeb(G_pos, 'w', 2);  % Plot the debundled edges with a specified width
    uistack(h_pos, 'top');  % Bring the node plot to the top

    % Customize the appearance of the positive edges
    set(he2_pos, 'FaceColor', "#bc0202");  % Set the edge color to a specific red shade

    % Add labels to nodes with custom styling
    labels = string(1:nModules);  % Create labels for each node (1 to nModules)
    for i = 1:nModules
        % Add labels to the top and bottom rows of nodes
        text(Xcorrds(i), 50, labels(i), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', 'white', 'FontSize', 8, ...
            'FontName', 'Helvetica', 'FontWeight', 'bold');
        text(Xcorrds(i), -50, labels(i), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', 'white', 'FontSize', 8, ...
            'FontName', 'Helvetica', 'FontWeight', 'bold');
    end

    % Final adjustments to the figure
    box off;  % Remove the box around the plot
    set(gca, 'Visible', 'off');  % Hide the axes
end
