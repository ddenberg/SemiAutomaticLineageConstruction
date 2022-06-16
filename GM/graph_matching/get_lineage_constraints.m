function lineage_constraints = get_lineage_constraints(lineage_corrections, node_IDs_1, node_IDs_2, current_frame, next_frame)

ind = ismember(lineage_corrections.frame_pairs, [current_frame, next_frame], 'rows');

node_IDs = lineage_corrections.node_pairs(ind,:);

node_map_1 = zeros(max(node_IDs_1), 1);
node_map_1(node_IDs_1) = 1:length(node_IDs_1);

node_map_2 = zeros(max(node_IDs_2), 1);
node_map_2(node_IDs_2) = 1:length(node_IDs_2);

lineage_constraints = [node_map_1(node_IDs(:,1)), node_map_2(node_IDs(:,2))];

end

