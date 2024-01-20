# An easy to use Bead detection tool for MATLAB

Do you have a lot of fluorescent beads and want to detect them accurately and quantify their time course?
Try BeadQuant.

The paper draft is under preparation. The user guide is `Bead 0.5.docx` in the top folder. 

## The BeadQuant GUI
The main interface.

![image](https://github.com/freemanwyz/BeadQuant/assets/3871340/1de6630e-c404-4ab7-b445-868379f1e236)

Detection results.

![image](https://github.com/freemanwyz/BeadQuant/assets/3871340/2804adf3-2735-4504-b0cd-be3f2af03c4a)

Quantify the time course.

![image](https://github.com/freemanwyz/BeadQuant/assets/3871340/0ec26c05-da80-4abe-b2c3-e04d98437ce8)


## Release notes
Version 0.5
	New user interface
	Support one or  two images
	Support one or two channels
	Select which channel is background and which is signal
	Support linear fitting
	Limit the range of τ and ‘a’ in curve fitting

Version 0.42  
	Option to choose radius range (in pixels)
	For each bead center, choose a circle with the highest mean value in channel 2. Extract the curve based on pixels in this curve as well as one circle outside it and one circle inside it.
	Bug fix: index for the control beads labeled in the output
 
Version 0.41  
	Bug fix: using unaligned background as F0
 
Version 0.4  
	- Select control beads after detection
	- Report curves, ΔF and τ for control beads
	- Save to Excel support platforms other than Windows
	- Include summary information in the Excel file
	- More output figure types
	- Bug fixes; Now tested on Matlab 2010b for better compatibility
 
Version 0.3  
	- Mean delta F from entire image (blue channel)
	- DeltaF and tau of the selected beads (blue channel)
	- Coordinates of location for selected beads
	- Time-lapse (in the blue channel) with the selected beads circled
	- Set Delta F/F or simply Delta F
 
Version 0.23  
	- Support LSM and TIFF
	- Can choose a time range in the GUI

	Algorithms updated to allow threshold setting for faster bead detection
	Set thresholds to accelerate the program for very dense image
	All the images are saved to the output folder as tiff format
