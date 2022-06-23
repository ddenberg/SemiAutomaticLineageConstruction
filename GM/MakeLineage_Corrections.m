

%Change path here to point to CPD2 folder
addpath(genpath('CPD2/core'));
addpath(genpath('CPD2/data'));
%path to graph_matching folder
addpath(genpath('graph_matching'))

% Name of registration output file
registration_filename = 'transforms.mat';
load(registration_filename);

% Name of graph matching output file
graph_output = 'graph_proposed.mat';
load(graph_output);

% Which pairs of frames to run over. Remember that the first frame is 0.
% If you would like to re-match certain frame pairs then set [frame_pairs] accordingly.
first_frame = 70;
final_frame = 72;
frame_pairs = [(first_frame:final_frame-1).', (first_frame+1:final_frame).'];

% also, check the alignment of this one with the time frame after
for ii = 1:size(frame_pairs, 1)
     
    % get pair of frames
    frame_pair = frame_pairs(ii,:);

    % Get index of registration struct
    registration_frame_pairs = cell2mat({registration.frame_pair}.');
    reg_ind = find(ismember(registration_frame_pairs, frame_pair, 'rows'));
    if isempty(reg_ind)
        error('Registration output not found for frame pair (%d, %d)', frame_pair(1), frame_pair(2));
    end

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
    
    
    nodes = [arrayfun(@(ind) sprintf('%03d_%03d', frame_pair(1), ind), uVal1, 'UniformOutput', false); 
             arrayfun(@(ind) sprintf('%03d_%03d', frame_pair(2), ind), uVal2, 'UniformOutput', false)];

    sample_graph = subgraph(G_lineage, nodes);
    
    node_pos_x = [centroids1(:,1); centroids2_transform(:,1)];
    node_pos_y = [centroids1(:,2); centroids2_transform(:,2)];
    node_pos_z = [centroids1(:,3); centroids2_transform(:,3)];

    node_colors = zeros(size(nodes, 1), 3);
    node_colors(1:size(centroids1, 1),1) = 1;
    node_colors(size(centroids1, 1)+1:end,3) = 1;

    % show point clouds registered (red is earlier time point)
    figure(1);
    clf;
    title('After registering Y to X.'); 
    scatter3(ptCloud1(:,1), ptCloud1(:,2), ptCloud1(:,3), '.r');
    hold on;
    scatter3(ptCloud2_transform(:,1), ptCloud2_transform(:,2), ptCloud2_transform(:,3), '.b');
    axis equal vis3d;
    title(sprintf('Pair (%d, %d)', frame_pair(1), frame_pair(2)));
    
    % visualization for checking if everything is correct - 3d plot of edges and nodes
    plot(sample_graph, 'XData', node_pos_x, 'YData', node_pos_y, 'ZData', node_pos_z, ...
        'EdgeColor', 'k', 'LineWidth', 2.0, 'Interpreter', 'none');

    figure(2);
    clf;
    plot(sample_graph, 'XData', node_pos_x, 'YData', node_pos_y, 'ZData', node_pos_z, ...
        'EdgeColor', 'k', 'LineWidth', 2.0,'NodeLabel',sample_graph.Nodes.Name, ...
        'NodeColor', node_colors, 'Interpreter', 'none');
    axis equal vis3d;
    title(sprintf('Pair (%d, %d)', frame_pair(1), frame_pair(2)));
        
    figure(3);
    clf;
    G_proposed = get_proposed_tree(G_lineage, frame_pair(2));
    plot(G_proposed,'layout','layered', 'Interpreter', 'none');
    title(sprintf('Pair (%d, %d)', frame_pair(1), frame_pair(2)));
    
    pause;
%     close all;
    
end

