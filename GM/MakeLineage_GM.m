

%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [G_lineage] is a MATLAB graph object containing 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Change path here to point to CPD2 folder
addpath(genpath('CPD2/core'));
addpath(genpath('CPD2/data'));
%path to graph_matching folder
addpath(genpath('graph_matching'))

% Name of registration output file
registration_filename = 'transforms.mat';
load(registration_filename);

% Graph output file
graph_output = 'graph_proposed.mat';

%Set graph input if you would like to restart from a previously constructed tree
graph_input = 'graph_proposed.mat';
if ~isempty(graph_input)
    load(graph_input);
else
    G_lineage = graph;
end

% flag to show plots after matching. You can view them again afterwards if you want.
show_plots = true;

% Which pairs of frames to run over. Remember that the first frame is 0.
% If you would like to re-match certain frame pairs then set [frame_pairs] accordingly.
first_frame = 70;
final_frame = 72;
frame_pairs = [(first_frame:final_frame-1).', (first_frame+1:final_frame).'];

%%%%%%%%%%%% IMPORTANT FLAG %%%%%%%%%%%%%%%%
%
% If you have made corrections to the lineage tree, set
% load_lineage_corrections to true to load them here. They will be used as
% constraints.
%

load_lineage_corrections = true;

lineage_corrections_filename = 'lineage_corrections.csv';
if load_lineage_corrections
    lineage_corrections = read_lineage_corrections(lineage_corrections_filename);
else
    lineage_corrections = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Graph matching parameter lambda. Lambda controls how much to weight cell
% volumes. If changing lambda does not produces better results then you will
% have to add manual corrections.
graph_match_lambda = 1;

% Lineage constraint. This controls how often cells are allowed to divide.
% If division_threshold = 20, then a daughter cell can only divide after 20
% frames after its mother divided.
division_time_threshold = 20;

% also, check the alignment of this one with the time frame after
for ii = 1:size(frame_pairs, 1)
    
    % get pair of frames
    frame_pair = frame_pairs(ii,:);

    fprintf('Beginning graph matching for pair (%d, %d)...\n', frame_pair(1), frame_pair(2));

    % Get index of registration struct
    registration_frame_pairs = cell2mat({registration.frame_pair}.');
    reg_ind = find(ismember(registration_frame_pairs, frame_pair, 'rows'));
    if isempty(reg_ind)
        error('Registration output not found for frame pair (%d, %d)', frame_pair(1), frame_pair(2));
    end
    
    % Get volumes
    volumes1 = registration(reg_ind).volumes1;
    volumes2 = registration(reg_ind).volumes2;

    % Compute average radius of a nucleus assuming they are spherical
    average_cell_volume = (sum(volumes1) + sum(volumes2)) / (length(volumes1) + length(volumes2));
    average_cell_radius = (average_cell_volume * 3 / (4 * pi))^(1/3);

    % Get transform between frames  
    Transform = registration(reg_ind).Transform;

    % Get centroids
    centroids1 = registration(reg_ind).centroids1;
    centroids2 = registration(reg_ind).centroids2;

    % Get point clouds
    ptCloud1 = registration(reg_ind).ptCloud1;
    ptCloud2 = registration(reg_ind).ptCloud2;

    % Get centroid labels
    uVal1 = registration(reg_ind).centroids1_ids;
    uVal2 = registration(reg_ind).centroids2_ids;
    
    % Transform ptCloud2 and centroids2
    ptCloud2_transform = cpd_transform(ptCloud2, Transform);
    centroids2_transform  = cpd_transform(centroids2, Transform);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sigma = 2 * average_cell_radius;
    A_old = squareform(pdist(centroids1));
    A_new = squareform(pdist(centroids2_transform));
    A_old = exp(-A_old.^2 / (2 * sigma^2));
    A_new = exp(-A_new.^2 / (2 * sigma^2));
    
    max_volume = max(max(volumes1), max(volumes2));
    B_old = volumes1 / max_volume;
    B_new = volumes2 / max_volume;

    Dist_new_old = pdist2(centroids1, centroids2_transform);
    X0 = exp(-Dist_new_old.^2 / (2 * sigma^2));
    
    if frame_pair(1) > 0
        G_lineage = get_proposed_tree(G_lineage, frame_pair(1));
        division_constraints = get_division_constraints(G_lineage, division_time_threshold, frame_pair(1));
    else
        division_constraints = [];
    end
    if ~isempty(lineage_corrections)
        lineage_constraints = get_lineage_constraints(lineage_corrections, uVal1, uVal2, frame_pair(1), frame_pair(2));
    else
        lineage_constraints = [];
    end

    if size(A_old, 1) ~= size(A_new, 1)
        X_min = ConstrainedMinimization_General(X0, A_old, A_new, B_old, B_new, graph_match_lambda, ...
            division_constraints, lineage_constraints);
    else
        X_min = FastPFP(X0, A_old, A_new, B_old, B_new, 0.5, graph_match_lambda, 1e-6, 1e-6, 1000, 100);
    end
    [~, match_ind] = max(X_min, [], 1);
    
    match_pairs = [uVal1(match_ind), uVal2];
    nodes = [arrayfun(@(ind) sprintf('%03d_%03d', frame_pair(1), ind), uVal1, 'UniformOutput', false); 
             arrayfun(@(ind) sprintf('%03d_%03d', frame_pair(2), ind), uVal2, 'UniformOutput', false)];
    
    edges = cell(size(match_pairs));
    edges(:,1) = arrayfun(@(ind) sprintf('%03d_%03d', frame_pair(1), ind), match_pairs(:,1), 'UniformOutput', false);
    edges(:,2) = arrayfun(@(ind) sprintf('%03d_%03d', frame_pair(2), ind), match_pairs(:,2), 'UniformOutput', false);
    edge_table = table(edges, 'VariableNames', {'EndNodes'});
    node_table_names = table(nodes, 'VariableNames', {'Name'});
    node_table_pos = table(nodes, [centroids1(:,1); centroids2_transform(:,1)], ...
                                  [centroids1(:,2); centroids2_transform(:,2)], ...
                                  [centroids1(:,3); centroids2_transform(:,3)], ...
                           'VariableNames', {'Name', 'xpos', 'ypos', 'zpos'});
    sample_graph = graph(edge_table, node_table_pos);
    
    G_lineage = addedge(G_lineage, edge_table);
    %Remove duplicate edges in lineage
    G_lineage = simplify(G_lineage);

    node_colors = zeros(size(nodes, 1), 3);
    node_colors(1:size(centroids1, 1),1) = 1;
    node_colors(size(centroids1, 1)+1:end,3) = 1;

    if show_plots
        % show point clouds registered (red is earlier time point)
        figure(1);
        clf;
        title('After registering Y to X.'); 
        scatter3(ptCloud1(:,1), ptCloud1(:,2), ptCloud1(:,3), '.r');
        hold on;
        scatter3(ptCloud2_transform(:,1), ptCloud2_transform(:,2), ptCloud2_transform(:,3), '.b');
        axis equal vis3d;
        
        % visualization for checking if everything is correct - 3d plot of edges and nodes
        plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, ...
            'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0, 'Interpreter', 'none');
        title(sprintf('Pair (%d, %d)', frame_pair(1), frame_pair(2)));
    
        figure(2);
        clf;
        plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, ...
            'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0,'NodeLabel',sample_graph.Nodes.Name, ...
            'NodeColor', node_colors, 'Interpreter', 'none');
        axis equal vis3d;
        title(sprintf('Pair (%d, %d)', frame_pair(1), frame_pair(2)));
            
        figure(3);
        clf;
        plot(G_lineage,'layout','layered', 'Interpreter', 'none');
        title(sprintf('Pair (%d, %d)', frame_pair(1), frame_pair(2)));
%         disp('time index');
%         disp(time_index);
        drawnow;
        pause;
    end

    fprintf(' Done!\n');

    % Save lineage tree
    save(graph_output, 'G_lineage');
end