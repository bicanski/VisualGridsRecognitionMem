This repo contains the code for my model of visual recogntion memory (Bicanski and Burgess 2019, Current Biology).

Bicanski A, Burgess N. - A computational model of recognition memory via grid cells. Current Biology, 2019, 29, 1–12. DOI: 10.1016/j.cub.2019.01.077

The code provided comes with no warranty etc. and you use it at your own risk. 
Furthermore it is licenced under the GPL. You are free to re-use and copy the 
code, but attribution must be correctly given.

I have tried focus on readbilty, not conciseness and/or efficiency. The code runs even on modest hardware as it 
is. I have done my best to use self-explanatory variable names and to annotate the code. 
If something remains cryptic I will happily assist. Please do not hesitate to contact me.

Andrej Bicanski, April 2019
andrej.bicanski@gmail.com

www.andrejbicanski.com



How to run the model:

1. Create grid cell firing rate maps with the file DOE_prep_makeGCs.m
You can also do this from the main script, see below

2. Set up some weights which will never change. This is done in the main script (DOE_main_SensoryPrediction.m), which has several flags at the very beginning, telling Matlab which bits of the file to execute.

A) setupDCs = 1; % set this to 1 to make grid cell/distance cell weights
B) setupGCs = 0; % set this to 1 to make grid cells as described above
C) getIMpts = 0; % set this to 1 to select the attention targets in images by hand
D) addIMpts = 0; % set this to 1 to add points to an image (e.g. distractors, see paper)
E) makeWTS  = 0; % set this to one learn the connections between sensory cells, feature label cells and grid cells. 

Weights created in A) are the same forever, i.e. they can be learned by an organism once, e.g. during development
Weights learned in E) are specific to your stimulus library, i.e. they correspond to one-shot learning when an agent encounters a new stimulus

3. Following this setup, the model's capability to recognize stimuli can be tested. Run DOE_main_SensoryPrediction.m 
again and set A-E above to zero and select the type of simulation from the list below:

cleanIMs = 1;   faceonly = 0;   % unperturbed stimuli
occluIMs = 0;   % stimuli with noise occlusions
realoIMs = 0;   % realistic occlusions (other stimuli)
famocIMs = 0;   % occlude with familiar stimuli
smallIMs = 0;   % downsized stimuli

