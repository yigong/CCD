ETCI: Electron Track Compton Imaging
====================================

### Purpose ###

Code for CCD image analysis, electron track reconstruction, electron track 
simulation, and other relevant stuff. bearing/ETCI can include any code you
have written which might be useful for someone else, now or in the future. At 
this point there is limited collaboration, in terms of people working together 
on the same code.

### Code practices ###

Let's try to move toward good code practices and style. Refer to Richard 
Johnson's ["The Elements of MATLAB Style"] [1], a copy of which Brian will 
soon have around his cubicle somewhere. Highlights from an earlier version are 
online at <http://www.datatool.com/downloads/matlab_style_guidelines.pdf>. That being 
said, Brian currently has lots of legacy codes which are stylistic 
catastrophes, which will still be in the repository for now because they are
important.

### Organization ###

The folders are intended to be used as follows.
- diffusion: For simulating detector response, usually from geant4 data.
- calibration\_and\_segmentation: Image processing from the raw \*.fit(s) data 
import, to energy calibration in each quadrant, to track segmentation.
- reconstruction: Algorithms for estimating the initial track direction, from 
a calibrated and segmented track image.
- visualization: Tools useful for viewing electron track images and simulated 
3D tracks.
- simulation: Archive of group member's working simulation code and simulation processing code.
- misc: generally useful codes, which might be required by other things

### Notes ###

If we have more languages of code than just Matlab, we should add subfolders 
for different languages, 'm' for matlab and 'py' for python or something.

[1]: <http://www.amazon.com/Elements-MATLAB-Style-Richard-Johnson/dp/0521732581/ref=sr_1_1?ie=UTF8&qid=1377218631&sr=8-1&keywords=matlab+style>
