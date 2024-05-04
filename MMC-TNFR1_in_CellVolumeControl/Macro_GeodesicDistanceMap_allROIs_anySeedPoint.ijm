/*
 * MACRO to create Geodesic Distance Map (GDM) between a small ROI (eg., a nucleus) or a point inside a bigger ROI (eg., a cell), for all ROIs in the ROI Manager
 * works on image stacks. ROIs should correspond to objects in every slice of stack
 * REQUIRED: MorphoLibJ plugin. Find 'IJPB-plugins' from Fiji>Update>Manage Update Sites
 * Cite the original paper in your work: https://doi.org/10.1093/bioinformatics/btw413
 * 
 * To use this macro...
 * 1. Add all Cell boundary ROIs in ROI Manager
 * 2. Add the corresponding Nucleus ROIs/centroids as point selections in ROI Manager IN THE SAME ORDER AS CELL ROIs
 * (if the main ROIs' centroid or center of mass is being used, no need for point 2.
 * 3. Click on the main image on which you want to perform the GDM analysis, run Macro
 */

//Check if the plugin in intalled
List.setCommands;
if (List.get("Geodesic Distance Map")=="") { exit("MorphoLibJ plugin not installed, install 'IJPB-plugins' from Fiji>Help>Update>Manage Update Sites"); }

//Get target image details
id=getTitle();
getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);

//Get user choice...
Choices=newArray("Nucleus ROI", "Geometric centroid of nucelus", "Manually designated points", "Geometric centroid of main ROI", "Intensity-weighted centroid (center of mass) of main ROI");
Dialog.create("Oeration Details...");
Dialog.addMessage("Choose seed for Godesic Distance Mapping (!!Options 1, 2, and 3 needs seed ROIs!!)...");
Dialog.addRadioButtonGroup("Select one", Choices, lengthOf(Choices), 1, "Geometric centroid of main ROI");
Dialog.show;
answer=Dialog.getRadioButton;
opMode=4;
for (i = 0; i < 5; i++) { if (answer==Choices[i]) { opMode=i+1; } }

nROIs=roiManager("count");
if (nROIs==0) { exit("No ROIs added!"); }

if (opMode<4)
{
	if ((nROIs%2)!=0) { exit("Number of Cell_ROIs is not equal to number of Nucleus_ROIs"); }
	nCellROIs=nROIs/2;
}
else { nCellROIs=nROIs; }

//Rename ROIs, get centroids
for (i = 0; i < nCellROIs; i++)
{
	roiManager("Select", i); roiManager("Rename", "MainROI_"+(i+1));
	if (opMode==1)
	{
		roiManager("Select", (i+nCellROIs)); roiManager("Rename", "NucleusROI_"+(i+1));
	}
	if (opMode==2)
	{
		roiManager("Select", (i+nCellROIs)); roiManager("Rename", "NucleusCentroidROI_"+(i+1));
		x=getValue("X raw"); y=getValue("Y raw");
		makePoint(x, y, "small yellow hybrid"); roiManager("add");
	}
	if (opMode==3)
	{
		roiManager("Select", (i+nCellROIs)); roiManager("Rename", "UserDefinedPoint_"+(i+1));
	}
	if (opMode==4)
	{
		roiManager("Select", i);
		x=getValue("X raw"); y=getValue("Y raw");
		makePoint(x, y, "small yellow hybrid"); roiManager("add");
	}
	if (opMode==5)
	{
		roiManager("Select", i);
		x=getValue("XM raw"); y=getValue("YM raw");
		makePoint(x, y, "small yellow hybrid"); roiManager("add");
	}
}

//Make a 32-bit blank image to copy and paste the GDMs in the next section
selectImage(id); run("Select None"); run("Duplicate...", "title=GDM duplicate");
selectImage("GDM"); run("Select None"); run("Multiply...", "value=0.0000000 stack"); run("32-bit");

//Make masks for each ROI, run GDM
for (i = 0; i < nCellROIs; i++)
{
	selectImage(id); run("Select None"); roiManager("Select", i); run("Create Mask");
	selectImage("Mask"); rename("MaskC");
	selectImage(id); run("Select None"); roiManager("Select", (i+nCellROIs)); run("Create Mask");
	selectImage("Mask"); rename("MaskN");
	run("Geodesic Distance Map", "marker=MaskN mask=MaskC distances=[Chessknight (5,7,11)] output=[32 bits] normalize");
	selectImage("MaskC-geoddist"); run("Select None"); roiManager("Select", i); run("Copy");
	selectImage("GDM"); run("Select None"); roiManager("Select", i); run("Paste");
	close("MaskC-geoddist"); close("MaskC"); close("MaskN");
}
selectImage("GDM"); run("Select None"); run("Duplicate...", "title=["+id+"_GDM] duplicate"); close("GDM");
selectImage(id+"_GDM"); run("Enhance Contrast", "saturated=0.35"); run("cool"); setSlice(1);