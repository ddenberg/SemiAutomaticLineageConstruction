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

In 'Preprocess_seg_errors.m' set the following variables:

1. [filename_seg_base] - path to the segmentation output in KLB format
2. [filename_raw_base] - path to the raw output in KLB format
