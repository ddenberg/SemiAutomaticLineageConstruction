# Graph Matching tool for Semi-Automatic Lineage Construction
MATLAB scripts:

**Preprocess_seg_errors.m**: Reads image segmentation output and identifies probable false positives regions and frames which may contain false negatives.

**Registration_Centroids.m**: Fast registration script which only uses cell centroids for the point cloud registration. ***NOTE: Right now only this output is compatible with the MakeLineage_GM.m script***.

**MakeLineage_GM.m**: Uses the registration output to create a lineage tree using sequential graph matching.

## Prerequisites
You must copy the MATLAB KLB wrapper function ***readKLBstack*** into this directory. You can find the KLB repository here: https://bitbucket.org/fernandoamat/keller-lab-block-filetype/src/master/.

If you are on Windows or Linux there are precompiled .mexw64 and .mexa64 files in the 'matlabWrapper' folder. If you are on MacOS you will need to compile the KLB library from source. I have included a separate MATLAB script, 'compileMex_mac.m', to compile the mex functions for MacOS after you have compiled the library. You will need to copy 'compileMex_mac.m' into the KLB directory in the folder 'matlabWrapper'.
