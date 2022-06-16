# SemiAutomaticLineageConstruction
Matlab code to register sequences of label images and then verify and fix lineage tracks as needed

Start by installing CPD for matlab https://sites.google.com/site/myronenko/research/cpd You need to fill out the form, download the stuff, if you don't have it already, install Xcode on your mac (this can take hours) 

### PrecomputeRegistrationTransforms.m
  This code will find the registration transforms between pairs of adjacent frames from the label images output by the 3D Stardist instance segmentation. You will need to change the first few lines of this code to point to your label images and to set the first and last frame that you want to register. The outputs will be put in the same folder as your label images.
  
### VisualizeSequence.m
  This code will create a video (.avi) file of the registered sequence using the outputs of the PrecomputeRegistrationTransforms. You can play this video using the open source VLC player.
