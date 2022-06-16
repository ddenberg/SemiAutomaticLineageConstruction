
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you don't have enough memory, processing power, or time to run the
% full registration, try running this script instead. It will do
% registration only on the cell centroids which still produces good
% results.
%
% It should take about 10-20 seconds per frame with the current parameters.
% This is runnable with 8GB of RAM but 16GB is recommended.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [store_registration] is a cell array of the best registration found over
% numTrials per pair of frames
%
% [store_sigma2] is a vector of the best sigma2 value found for each
% registered pair
%
% [store_centroids] is a cell array of the unregistered centroids between each pair of frames
%
% [store_point_clouds] is a cell array of the unregistered point clouds between each pair of frames
%
% [store_volumes] is a cell array of the volumes for each nucleus between each pair of frames
%
% Remember that index 1 of these arrays corresponds to the pair of frames
% [0, 1] (more generally index n -> frame pair [n-1, n]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set numThreads to the number of cores in your computer. If your processor
% supports hyperthreading/multithreading then set it to 2 x [number of cores]
numThreads = 4;

addpath(genpath('CPD2/core'));
addpath(genpath('CPD2/data'));

% What is the prefix for the embryo names?
% filename_seg_base = 'E:/Posfai_Lab/rpky/220309_out/st0/klb/klbOut_Cam_Long_%05d.lux.label.klb';
filename_seg_base = '/media/david/Seagate_Exp/Posfai_Lab/rpky/220309_out/st0/klb/klbOut_Cam_Long_%05d.lux.label.klb';

% Name of output file
Registration_filename = 'transforms.mat';

% Set this to false if you do not want to use the preprocessing output
use_preprocess_false_positives = true;

% Which image frames to run over. Remember that the first frame is 0
valid_time_indices = 0:100;

% Voxel size before making isotropic
pixel_size_xy_um = 0.208; % um
pixel_size_z_um = 2.0; % um
% Voxel size after making isotropic
xyz_res = 0.8320;
% Volume of isotropic voxel
voxel_vol = xyz_res^3;

% Threshold to accept registration
sigma2_threshold = 5;

% Name of preprocessing output
Preprocess_filename = 'preprocess.mat';
if use_preprocess_false_positives
    preprocess = load(Preprocess_filename);
else
    preprocess = struct('store_false_positives_guess', cell(length(valid_time_indices), 1));
end

% Option to downsample point clouds in the saved output. Will reduce storage space
downsample_fraction = 1/10;

% Initialize empty graph and cell array for storing registration
store_centroids = cell(length(valid_time_indices) - 1, 2);
store_registration = cell((length(valid_time_indices)-1), 1);
store_sigma2 = zeros((length(valid_time_indices)-1), 1);
store_point_clouds = cell(length(valid_time_indices) - 1, 2);
store_volumes = cell(length(valid_time_indices) - 1, 2);

% numTrials controls how many random initializations to do 
numTrials = 1e3;

tic;
for time_index_index = 1:length(valid_time_indices)-1 
     
    fprintf('Beginning Registration index %d...', time_index_index);
    
    % time_index is the frame id of the current frame
    time_index = valid_time_indices(time_index_index);
    
    % time_index_plus_1 is the frame id of the following frame
    time_index_plus_1 = valid_time_indices(time_index_index+1);
    
    % read in segmented images
    filename_seg = sprintf(filename_seg_base, time_index);
    seg1 = readKLBstack(filename_seg, numThreads);
    
    filename_seg = sprintf(filename_seg_base, time_index_plus_1);
    seg2 = readKLBstack(filename_seg, numThreads);
    
    % Exclude regions in segmented image whose mean intensities are within 2 stds of the background
    if use_preprocess_false_positives
        seg1(ismember(seg1, preprocess.store_exclude{time_index_index})) = 0;
        seg2(ismember(seg2, preprocess.store_exclude{time_index_index+1})) = 0;
    end
    
    % Rescale image 
    resXY = 0.208;
    resZ = 2.0;
    reduceRatio = 1/4;
    seg1 = isotropicSample_nearest(seg1, resXY, resZ, reduceRatio);
    seg2 = isotropicSample_nearest(seg2, resXY, resZ, reduceRatio);
   
    % Find non-zero indices of image
    [Y, XZ, Val1] = find(seg1);
    [X, Z] = ind2sub([size(seg1, 2), size(seg1, 3)], XZ);
    
    ptCloud1 = [X, Y, Z];
    ptCloud1 = ptCloud1 - mean(ptCloud1, 1);
    
    %Compute centroids
    uVal1 = unique(Val1);
    centroids1 = zeros(length(uVal1), 3);
    for ii = 1:length(uVal1)
        centroids1(ii,:) = mean(ptCloud1(Val1 == uVal1(ii),:), 1);
    end
    
    [Y, XZ, Val2] = find(seg2);
    [X, Z] = ind2sub([size(seg2, 2), size(seg2, 3)], XZ);
    
    ptCloud2 = [X, Y, Z];
    ptCloud2 = ptCloud2 - mean(ptCloud2, 1);
    
    %Compute centroids
    uVal2 = unique(Val2);
    centroids2 = zeros(length(uVal2), 3);
    for ii = 1:length(uVal2)
        centroids2(ii,:) = mean(ptCloud2(Val2 == uVal2(ii),:), 1);
    end

    % Compute some stats about the volumes of each cell to make an estimate of the cell radius
    stats1 = regionprops3(seg1, 'Volume');
    volumes1 = stats1.Volume(stats1.Volume > 0);
    stats2 = regionprops3(seg2, 'Volume');
    volumes2 = stats2.Volume(stats2.Volume > 0);

    % optionally downsample point clouds
    if downsample_fraction < 1
        p1 = randperm(size(ptCloud1, 1), round(size(ptCloud1, 1) * downsample_fraction));
        p2 = randperm(size(ptCloud2, 1), round(size(ptCloud2, 1) * downsample_fraction));
        
        find1 = find(seg1);
        find2 = find(seg2);

        find1 = find1(p1);
        find2 = find2(p2);

        Val1 = Val1(p1);
        Val2 = Val2(p2);

        ptCloud1 = ptCloud1(p1,:);
        ptCloud2 = ptCloud2(p2,:);
    end


    moving = pointCloud(centroids2);
    fixed = pointCloud(centroids1);

    step = 1;
    sigma2 = Inf;
    sigma2_trials = zeros(numTrials,1);
    transforms_trials = cell(numTrials, 1);
    while ((sigma2 > sigma2_threshold) && (step <= numTrials))    

        if step == 1
            R = eye(3);
        else
            R = orth(randn(3));
        end

        tform = rigid3d(R, [0, 0, 0]);
        
        moving_temp = pctransform(moving,tform);

        % Set the options
        opt.method='rigid'; % use rigid registration
        opt.viz=0;          % show every iteration
        opt.outliers=0;     % do not assume any noise

        opt.normalize=0;    % normalize to unit variance and zero mean before registering (default)
        opt.scale=0;        % estimate global scaling too (default)
        opt.rot=1;          % estimate strictly rotational matrix (default)
        opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)

        opt.max_it=200;     % max number of iterations
        opt.tol=1e-5;       % tolerance
        opt.ftg = 1; % make faster


        % registering Y to X
        [Transform, ~, sigma2] = cpd_register(fixed.Location, moving_temp.Location, opt);

        Transform.R = (R * Transform.R')';
        
        store_registration{time_index_index, 1} = Transform;
        store_sigma2(time_index_index) = sigma2;
        
        transforms_trials{step, 1} = Transform;
        sigma2_trials(step) = sigma2;
        step = step + 1;
    end

    if step > numTrials
%         fprintf(' Did not find transformation with sigma2 < %f', sigma2_threshold);
        % get the best one we found
        [min_sigma2, min_ind] = min(sigma2_trials);
        Transform = transforms_trials{min_ind,1};
        store_registration{time_index_index, 1} = Transform;
        store_sigma2(time_index_index) = min_sigma2;
    end

    store_centroids{time_index_index,1} = centroids1;
    store_centroids{time_index_index,2} = centroids2;

    store_point_clouds{time_index_index,1} = ptCloud1;
    store_point_clouds{time_index_index,2} = ptCloud2;

    store_volumes{time_index_index,1} = volumes1;
    store_volumes{time_index_index,2} = volumes2;

    fprintf(' Best Sigma2: %f, Done!\n', store_sigma2(time_index_index));
end
toc;

% Save output
save(Registration_filename, 'store_registration', 'store_sigma2', 'store_centroids', ...
    'store_point_clouds', 'store_volumes');
