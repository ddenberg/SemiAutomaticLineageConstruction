# Graph Matching tool for Semi-Automatic Lineage Construction
MATLAB scripts:

**Preprocess_seg_errors.m**: Reads image segmentation output and identifies probable false positives regions and frames which may contain false negatives.

**Registration_Centroids.m**: Fast registration script which only uses cell centroids for the point cloud registration. ***NOTE: Right now only this output is compatible with the MakeLineage_GM.m script***.

**MakeLineage_GM.m**: Uses the registration output to create a lineage tree using sequential graph matching.

**MakeLineage_Corrections.m**: Visualize the graph matching output to make corrections where necessary.

## Prerequisites
You must copy the MATLAB KLB wrapper function ***readKLBstack*** into this directory. You can find the KLB repository here: https://bitbucket.org/fernandoamat/keller-lab-block-filetype/src/master/.

If you are on Windows or Linux there are precompiled .mexw64 and .mexa64 files in the 'matlabWrapper' folder. If you are on MacOS you will need to compile the KLB library from source. I have included a separate MATLAB script, 'compileMex_mac.m', to compile the mex functions for MacOS after you have compiled the library. You will need to copy 'compileMex_mac.m' into the KLB directory in the folder 'matlabWrapper'.

### Compiling The KLB Library on Mac or Linux
Ideally you won't have to do this, but if the mex files do not work for whatever reason, or you are on a Mac then follow these steps:

1. You will need to have Cmake and a C++ compiler installed on your computer.
2. Navigate to the KLB main directory.
3. Open a terminal and enter the following
```
mkdir build
cd build
cmake ..
make
```

If you have gotten no errors, then the KLB library should have been compiled successfully. You can now copy 'compileMex_mac.m' to the matlabWrapper folder and compile the 'readKLBstack' mex file.

## Preprocessing Instructions

In 'Preprocess_seg_errors.m' set the following parameters:

1. [filename_seg_base] - path to the segmentation output in KLB format
2. [filename_raw_base] - path to the raw output in KLB format
3. [output_name] - name of .mat file where the preprocessing output will be saved.
4. [final_frame] - last frame ID of the image to process
5. [background_std_threshold] - how many standard deviations of the background noise to flag false positives
6. [do_false_negatives_filter] - flag to filter for false negatives (uses more memory)
7. [volume_threshold] - volume threshold for false negatives (does not matter if [do_false_negatives_filter] is false)
8. [cell_std_threshold] - how many standard deviations of the nuclei signal to flag for false negatives

This script will output two arrays which will be stored in [output_name]:

1. [store_false_positives_guess] - cell array which stores the IDs of regions identified to be false positives for each frame
2. [store_false_negatives_guess] - logical array which is True if there may be a false negative at the corresponding index
3. [store_numcells] - number of cells per frame with probable false positives excluded

Remember that index 1 of  these arrays corresponds to frame 0 in the image files. (Index n -> Frame n-1)

### Formatting [filename_seg_base] and [filename_raw_base]

Make sure to use formatting operators to specify the ID of each frame. If you include more than one formatting operator know that they will all be filled with the same value.

Example: 
```
filename_raw_base = '[LOCATION]/stack_0_channel_0_obj_left/out/folder_Cam_Long_%05d.lux/klbOut_Cam_Long_%05d.lux.klb';
```
In this case I know that the raw images will look like 'stack_0_channel_0_obj_left/out/folder_Cam_Long_00023.lux/klbOut_Cam_Long_00023.lux.klb' with up to 5 leading zeros where each ID appears.

### Correcting Segmentation Errors

From the list of probable false positives and false negative frames you can now use Abhishek's modified AnnotatorJ tool to correct regions of interest. If you find that the false positives were correctly identified (or trust that they were correctly identified) then you can use a flag in Registration_Centroids.m to automatically ignore those regions.
