
%%%%%%%%%%%%%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script may use a lot of memory as it involves reading the raw image
% file. ~8GB of memory is recommended.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [preprocess] array of structs which contains:
%       -frame_id: ID of the frame
%       -false_positive_ids_guess: list of region ids identified to be false positives
%       -false_negative_guess: True if a false negative is found
%       -num_cells: guess of number of cells per frame (after excluding the false positive guesses)
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
preprocess_output = 'preprocess.mat';
if isfile(preprocess_output)
    load(preprocess_output);
else
    preprocess = [];
end

% Which image frames to run over. Remember that the first frame is 0
final_frame = 100;
valid_frames = 0:final_frame;

% How many standard deviations within background noise to flag false positives
background_std_threshold = 2;

% Do false negatives (set this to false to preserve memory)
do_false_negatives_filter = false;

%Volume threshold for false negatives. If you are overestimating the number
%of false negatives, try increasing this value.
volume_threshold = 3000;

% How many standard deviations within nuclei signal to flag for false negative
cell_std_threshold = 2;

store_false_positives_guess = cell(length(valid_frames), 1);
store_numcells = zeros(length(valid_frames), 1);
store_false_negatives_guess = false(length(valid_frames), 1);

tic;
for ii = 1:length(valid_frames) 

    % get frame id
    frame_id = valid_frames(ii);

    fprintf('Beginning Preprocessing Frame ID: %d...', frame_id);
    
    % read in segmented images
    filename_seg_base_nspec = count(filename_seg_base, '%');
    filename_seg = sprintf(filename_seg_base, frame_id * ones(1, filename_seg_base_nspec));
    seg = readKLBstack(filename_seg, numThreads);
    
    % read in raw images
    filename_raw_base_nspec = count(filename_raw_base, '%');
    filename_raw = sprintf(filename_raw_base, frame_id * ones(1, filename_raw_base_nspec));
    raw = readKLBstack(filename_raw, numThreads);
    
    % Exclude regions in segmented image whose mean intensities are within 2 stds of the background
    foreground_ind = find(seg > 0);
    num_bg = numel(seg) - length(foreground_ind);
    s_all = sum(raw, 'all');
    s_foreground = sum(raw(foreground_ind), 'all');
    background_mean = (s_all - s_foreground) / num_bg;

    ss_all = double(sum((single(raw) - background_mean).^2, 'all'));
    ss_foreground = double(sum((single(raw(foreground_ind)) - background_mean).^2, 'all'));
    ss_bg = ss_all - ss_foreground;
    background_std = sqrt(ss_bg / (num_bg - 1));

    stats = regionprops3(seg, raw, {'MeanIntensity', 'Volume'});
    exclude_logical = abs(stats.MeanIntensity - background_mean) < background_std_threshold * background_std;
    exclude_id = find(exclude_logical);
    
    nan_ind = isnan(stats.MeanIntensity);
    numcells = length(stats.MeanIntensity(~exclude_logical & ~nan_ind));
    cell_intensity_mean = mean(stats.MeanIntensity(~exclude_logical & ~nan_ind));
    cell_intensity_std = std(stats.MeanIntensity(~exclude_logical & ~nan_ind));
    
    frame_false_negative = false;
    if do_false_negatives_filter
        raw_subtract = raw;
        raw_subtract(imdilate(seg > 0, strel('sphere', 10))) = background_mean;
        BW = raw_subtract > (background_mean + 2 * background_std);
        BW = bwareaopen(BW, 100, 26);
        
        stats_BW = regionprops3(BW, raw_subtract, {'Volume', 'MeanIntensity'});
        if any(stats_BW.Volume > volume_threshold & ...
               abs(stats_BW.MeanIntensity - cell_intensity_mean) < cell_std_threshold * cell_intensity_std)
            frame_false_negative = true;
        end
    end

    % Update output struct
    if ~isempty(preprocess)
        stored_frame_ids = [preprocess.frame_id].';
        ind = find(ismember(stored_frame_ids, frame_id, 'rows'));

        if isempty(ind)
            ind = size(preprocess, 1) + 1;
        end

        preprocess(ind,1) = struct('frame_id', frame_id, ...
                                   'false_positive_ids_guess', exclude_id, ...
                                   'num_cells', numcells, ...
                                   'false_negative_guess', frame_false_negative);
    else
        preprocess = struct('frame_id', frame_id, ...
                            'false_positive_ids_guess', exclude_id, ...
                            'num_cells', numcells, ...
                            'false_negative_guess', frame_false_negative);
    end 
    
    fprintf(' Done!\n');
end
toc;

% Loop through frames and identify if the number of cells goes down between frames
stored_frame_ids = [preprocess.frame_id].';
for ii = 1:size(preprocess, 1)
    frame_id = preprocess(ii).frame_id;
    if frame_id > 0
        frame_id_prev = frame_id - 1;
        prev_ind = find(stored_frame_ids == frame_id_prev);
	    if ~isempty(prev_ind)
		    if preprocess(ii).num_cells < preprocess(prev_ind).num_cells
                preprocess(ii).false_negative_guess = true;
            end
		end
    end
end

save(preprocess_output, 'preprocess');