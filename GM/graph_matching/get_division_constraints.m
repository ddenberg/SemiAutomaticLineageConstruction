function division_constraints = get_division_constraints(G_lineage, time_threshold, current_frame)
% time_threshold is the number of steps backward in the tree to look for a division
% the return vector division_constraints contain the nodes which come from a division within the time_threshold

%Get nodes in G_lineage
nodes = G_lineage.Nodes.Name;
nodes = cellfun(@(str) split(str, '_'), nodes, 'UniformOutput', false);
nodes = horzcat(nodes{:}).';
nodes = cellfun(@(str) str2double(str), nodes, 'UniformOutput', false);
nodes = cell2mat(nodes);
%first column of nodes is the frame index and the second column is the node id in that frame

%get terminal nodes at the current frame
terminal_nodes = nodes(:,1) == current_frame;
terminal_nodes_ind = find(terminal_nodes);

terminal_node_names = arrayfun(@(ind) sprintf('%03d_%03d', nodes(ind,1), nodes(ind,2)), terminal_nodes_ind, 'UniformOutput', false);

% Make sure we are working with the proposed tree
G_lineage = get_proposed_tree(G_lineage, current_frame);

%compute the distances between all the terminal nodes
D = distances(G_lineage, terminal_node_names, terminal_node_names);

%Pairs of nodes whose distance is > 2*time_threshold have a shared parent longer than time_threshold away
D_threshold = D < 2 * time_threshold;

ind = find(triu(D_threshold, 1));
[I,J] = ind2sub(size(D), ind);

division_constraints = unique([I;J]);

end

