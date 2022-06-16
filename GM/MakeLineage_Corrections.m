

%Change path here to point to CPD2 folder
addpath(genpath('CPD2/core'));
addpath(genpath('CPD2/data'));
%path to graph_matching folder
addpath(genpath('graph_matching'))

% Name of registration output file
registration_filename = 'transforms.mat';
registration_data = load(registration_filename);

% Name of graph matching output file
graph_output = 'graph_proposed.mat';
load(graph_output);

% Which image frames to run over. Remember that the first frame is 0
final_frame = 40;
valid_time_indices = 1:final_frame;

% also, check the alignment of this one with the time frame after
for time_index_index = 1:(length(valid_time_indices)-1)
     
    % store this time index
    time_index = valid_time_indices(time_index_index);
    
    % store next in series
    time_index_plus_1 = valid_time_indices(time_index_index+1);

    % Get transform between frames  
    Transform = registration_data.store_registration{time_index_plus_1,1};

    % Get centroids
    centroids1 = registration_data.store_centroids{time_index_plus_1,1};
    centroids2 = registration_data.store_centroids{time_index_plus_1,2};

    % Get point clouds
    ptCloud1 = registration_data.store_point_clouds{time_index_plus_1,1};
    ptCloud2 = registration_data.store_point_clouds{time_index_plus_1,2};

    % Get centroid labels
    uVal1 = registration_data.store_centroid_ids{time_index_plus_1,1};
    uVal2 = registration_data.store_centroid_ids{time_index_plus_1,2};
    
    % Transform ptCloud2 and centroids2
    ptCloud2_transform = cpd_transform(ptCloud2, Transform);
    centroids2_transform  = cpd_transform(centroids2, Transform);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    nodes = [arrayfun(@(ind) sprintf('%03d_%03d', time_index, ind), uVal1, 'UniformOutput', false); 
             arrayfun(@(ind) sprintf('%03d_%03d', time_index_plus_1, ind), uVal2, 'UniformOutput', false)];

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
    
    % visualization for checking if everything is correct - 3d plot of edges and nodes
    plot(sample_graph, 'XData', node_pos_x, 'YData', node_pos_y, 'ZData', node_pos_z, ...
        'EdgeColor', 'k', 'LineWidth', 2.0, 'Interpreter', 'none');

    figure(2);
    clf;
    plot(sample_graph, 'XData', node_pos_x, 'YData', node_pos_y, 'ZData', node_pos_z, ...
        'EdgeColor', 'k', 'LineWidth', 2.0,'NodeLabel',sample_graph.Nodes.Name, ...
        'NodeColor', node_colors, 'Interpreter', 'none');
    axis equal vis3d;
        
    figure(3);
    clf;
    G_proposed = get_proposed_tree(G_lineage, time_index_plus_1);
    plot(G_proposed,'layout','layered', 'Interpreter', 'none');
    disp('time index');
    disp(time_index);
    
    pause;
%     close all;
    
end

