
%%%%%%%%%%%%%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script may use a lot of memory as it involves reading the raw image
% file. ~16GB of memory is recommended.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [store_false_positives_guess] is a cell array which contains the IDs of
% labels which this script has identified to be false positives. Note: the
% ID corresponds to the ID of the label in the klb file (i.e. that is the
% number you will see in imageJ).
%
% [store_false_negatives] is a logical array which is True if there may be
% a false negative in the corresponding frame.
%
% Remember that index 1 of these arrays corresponds to frame 0 in the image
% files. (index n -> frame n-1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set numThreads to the number of cores in your computer. If your processor
% supports hyperthreading/multithreading then set it to 2 x [number of cores]
numThreads = 4;

% What is the prefix for the embryo names?
% filename_seg_base = 'E:/Posfai_Lab/rpky/220309_out/st0/klb/klbOut_Cam_Long_%05d.lux.label.klb';
% filename_raw_base = 'E:/Posfai_Lab/rpky/220309/stack_0_channel_0_obj_left/out/folder_Cam_Long_%05d.lux/klbOut_Cam_Long_%05d.lux.klb';
filename_seg_base = '/media/david/Seagate_Exp/Posfai_Lab/rpky/220309_out/st0/klb/klbOut_Cam_Long_%05d.lux.label.klb';
filename_raw_base = '/media/david/Seagate_Exp/Posfai_Lab/rpky/220309/stack_0_channel_0_obj_left/out/folder_Cam_Long_%05d.lux/klbOut_Cam_Long_%05d.lux.klb';

% Change name or destination of the output
output_name = 'preprocess.mat';

pixel_size_xy_um = 0.208; % um
pixel_size_z_um = 2.0; % um
% Voxel size after making isotropic
xyz_res = 0.8320;
% Volume of isotropic voxel
voxel_vol = xyz_res^3;

%Volume threshold for false negatives. If you are overestimating the number
%of false negatives, try increasing this value.
volume_threshold = 3000;

% Which image frames to run over. Remember that the first frame is 0
valid_time_indices = 0:100;

% Do false negatives (set this to false to preserve memory)
do_false_negatives_filter = false;

store_false_positives_guess = cell(length(valid_time_indices), 1);
store_numcells = zeros(length(valid_time_indices), 1);
store_false_negatives = false(length(valid_time_indices), 1);

for ii = 1:length(valid_time_indices)     
    fprintf('Beginning Preprocessing index %d...', ii);
    
    % store this time index
    time_index = valid_time_indices(ii);
    
    % read in segmented images
    filename_seg = sprintf(filename_seg_base, time_index);
    seg = readKLBstack(filename_seg, numThreads);
    
    % read in raw images
    filename_raw = sprintf(filename_raw_base, time_index, time_index);
    raw = readKLBstack(filename_raw, numThreads);
    
    % Exclude regions in segmented image whose mean intensities are within 2 stds of the background
    background_ind = seg == 0;
    background_mean = mean(raw(background_ind));
    background_std = std(single(raw(background_ind)));
    stats = regionprops3(seg, raw, {'MeanIntensity', 'Volume'});
    exclude_logical = abs(stats.MeanIntensity - background_mean) < 2 * background_std;
    exclude_id = find(exclude_logical);
    
    numcells = length(stats.MeanIntensity(~exclude_logical));
    cell_intensity_mean = mean(stats.MeanIntensity(~exclude_logical));
    cell_intensity_std = std(stats.MeanIntensity(~exclude_logical));
    
    if do_false_negatives_filter
        raw_subtract = raw;
        raw_subtract(imdilate(seg > 0, strel('sphere', 10))) = background_mean;
        BW = raw_subtract > (background_mean + 2 * background_std);
        BW = bwareaopen(BW, 100, 26);
        
        stats_BW = regionprops3(BW, raw_subtract, {'Volume', 'MeanIntensity'});
        if any(stats_BW.Volume > volume_threshold & ...
               abs(stats_BW.MeanIntensity - cell_intensity_mean) < 2 * cell_intensity_std)
            store_false_negatives(ii) = true;
        end
    end

    if ii > 1
        if numcells < store_numcells(ii-1)
            store_false_negatives(ii) = true;
        end
    end
    
    store_false_positives_guess{ii} = exclude_id;
    store_numcells(ii) = numcells;
    
    fprintf(' Done!\n');
end

save(output_name, 'store_false_positives_guess', 'store_false_negatives_guess');