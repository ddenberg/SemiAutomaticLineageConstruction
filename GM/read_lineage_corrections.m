function lineage_corrections = read_lineage_corrections(filename)

raw = readmatrix(filename, 'Delimiter', ',', 'OutputType', 'char');

num_corrections = size(raw, 1);

frame_pairs = zeros(num_corrections, 2);
node_pairs = zeros(num_corrections, 2);
for ii = 1:num_corrections
    temp = split(raw{ii,1}, '_');
    frame_pairs(ii,1) = str2double(temp{1});
    node_pairs(ii,1) = str2double(temp{2});

    temp = split(raw{ii,2}, '_');
    frame_pairs(ii,2) = str2double(temp{1});
    node_pairs(ii,2) = str2double(temp{2});

end

lineage_corrections.frame_pairs = frame_pairs;
lineage_corrections.node_pairs = node_pairs;
end

