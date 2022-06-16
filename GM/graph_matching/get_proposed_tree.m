function G_proposed = get_proposed_tree(G_lineage, current_frame)

%Get nodes in G_lineage
nodes = G_lineage.Nodes.Name;
nodes = cellfun(@(str) split(str, '_'), nodes, 'UniformOutput', false);
nodes = horzcat(nodes{:}).';
nodes = cellfun(@(str) str2double(str), nodes, 'UniformOutput', false);
nodes = cell2mat(nodes);
%first column of nodes is the frame index and the second column is the node id in that frame

% Remove nodes in graph that are after the current_frame
remove_nodes = nodes(:,1) > current_frame;
remove_nodes_ind = find(remove_nodes);
remove_node_names = arrayfun(@(ind) sprintf('%03d_%03d', nodes(ind,1), nodes(ind,2)), remove_nodes_ind, 'UniformOutput', false);

G_proposed = rmnode(G_lineage, remove_node_names);
end

