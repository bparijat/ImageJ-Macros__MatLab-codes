/*
 * for FIJI
 * Macro to open BEAD image stack and run Anisotropy processing on all images in folder
 * for a set of N number of images/stacks, if there are N BEAD images, put all BEAD images in a stack SEQUENTIALLY
 * THERE SHOULD BE NO OTHER IMAGES IN THE FOLDER TO BE ANALYSED
 * name BEAD image stack as aaa000BEADS.tif and save images in the same folder where images to be analysed are kept
 * MAKE SURE BEAD IMAGE STACK IS AT TOP WHEN IMAGES IN FOLDER ARE NAMED IN ASCENDING ORDER
 * images must be arranged sequentially as bead images, so name them accordingly
 * run macro, when prompted, browse to BEAD image stack from folder, select it
 * when prompted, select folder where images to be analysed are kept
 * For proper values of G-factor, an image whose name starts with the word "Gfactor" is needed to calculate correction factor. The image must have zero anisotropy. Or else a single number for Gfactor can be used can be used
 * Have the Gfactor image open
 * PROCESSING DETAILS
 * LEFT half of image=PERPENDICULAR channel, RIGHT half of image=PARALLEL channel
 * original code has been written for 2048x2048 image
 * 
 * Plugins required: 1. Descriptor-based registration (2d/3d)  2. Linear Stack Alignment With SIFT
 * 
 * image of sub-pixel beads is registered using Descriptor-based registration (2d/3d)
 * using the same transformation matrix for beads, cell images are registered and then anisotropy is calculated
 * SIFT Registration is a backup
 * 
 */

//FOR USER INPUT
//savedfilepath="C:\\Users\\IACS\\Desktop\\AnisotropyNew\\";
//savedfilepath="C:\\Users\\Dell\\Desktop\\Anisotropy\\";
//savedfilepath="C:\\Users\\Lord Parijat\\Desktop\\Anisotropy\\";

dir = getDirectory("Browse to Folder where bead image stack and images to be analysed are kept");
list = getFileList(dir);
nf=list.length;
savedfilepath=dir;


beadsz=5; //approximate size of beads in pixels
backPerpendicular=100;
minAni=0.15; maxAni=0.25;  //scale of Anisotropy values that should give the best contrast
gFactorUNPolarised=0.67; //for bead image, if taken in same polarisation channel

//give file names of bead and Gfactor images
GfactorImage="None";
nOpenImgs=nImages;
for (imgID1=0; imgID1<nOpenImgs; imgID1++)
{
selectImage(imgID1+1); imageFileName=getTitle(); ltt=lengthOf(imageFileName);
if (ltt>=7) { gstr=substring(imageFileName, 0, 7); if (gstr=="Gfactor") { GfactorImage=imageFileName; } }
}

xshift=1024; //round(cellWid/2); //the perpendicular channel being 'xshift' pixels away horizontally
yshift=0; //if there is a vertical shift between the 2 channels
xCutOff=30; yCutOff=xCutOff; //pixels to be cropped from X and Y axis (to avoid discontinuation in images)

//marker flag, put pflag=1 to save parallel and perpendicular channels separately
//else put pflag=0 to not save them
pflag=0;

//marker flag, put StRegflag=1 to use StackReg/SIFT plugin for image registration
//else put StRegflag=0 to not use it
StRegflag=0;

//marker flag, put maskcheck=1 to use a mask to isolate cell from background
//else put maskcheck=0 to not use it
maskcheck=0;

Gfactor=0.67;
Dialog.create("Operation Details");
if (GfactorImage=="None")
{
	Dialog.addMessage("No G-factor image is open");
	Dialog.addNumber("Use this G-factor value:", Gfactor); gfm=0;
	Dialog.addMessage("Reference:\nGreen Fluorophore-Objective-10x=0.680, 20x=0.683, 40x=0.660, 63x=0.670\nor Red Fluorophore-Objective-40x=0.664, 63x=0.700");
}
else { Dialog.addMessage("Calculate G-factor with this Image="+GfactorImage); gfm=1; }
Dialog.addChoice("Transformation Model for Image Registration:", newArray("Homography", "Affine", "Rigid"));
Dialog.addNumber("Average Background (left channel)", backPerpendicular);
Dialog.addNumber("Size of beads [pixels] (for image registration)", beadsz);
Dialog.addNumber("Horizontal Distance between channels (pixels))", xshift);
Dialog.addNumber("Vertical Distance between channels (pixels)", yshift);
Dialog.addMessage("To avoid image discontinuity in final image,");
Dialog.addNumber("Pixels to be cropped from X and Y axis", xCutOff);
Dialog.addCheckbox("Use Intensity Mask (AutoThreshold)", false);
Dialog.addCheckbox("Additionally use Linear Stack Alignment With SIFT", false);
Dialog.addCheckbox("Save registered images separately", false);
Dialog.addNumber("Minimum Anisotrpy displayed:", minAni);
Dialog.addNumber("Maximum Anisotrpy displayed:", maxAni);
Dialog.show();

if (GfactorImage=="None") { Gfactor=Dialog.getNumber(); }
TransformModel=Dialog.getChoice();
backPerpendicular=Dialog.getNumber();
beadsz=Dialog.getNumber();
xshift=Dialog.getNumber();
yshift=Dialog.getNumber();
xCutOff=Dialog.getNumber(); yCutOff=xCutOff;
chkbx1=Dialog.getCheckbox();
chkbx2=Dialog.getCheckbox();
chkbx3=Dialog.getCheckbox();
if (chkbx1==true)  { maskcheck=1; }
if (chkbx1==false)  { maskcheck=0; }
if (chkbx2==true)  { StRegflag=1; }
if (chkbx2==false)  { StRegflag=0; }
if (chkbx3==true)  { pflag=1; }
if (chkbx3==false)  { pflag=0; }
minAni=Dialog.getNumber();
maxAni=Dialog.getNumber();

backParallel=0.8*backPerpendicular;


//identify the bead image stack in the folder
beadImgIndex=99999999;
for (i=0; i<nf; i++)
{
	open(dir+list[i]);
	beadImgMM=getTitle();
	if (beadImgMM=="aaa000BEADS.tif") { beadImgIndex=i; break; }
	else { selectWindow(beadImgMM); close(); continue; }
}

if (beadImgIndex==99999999)
{
	Dialog.create("WARNING");
	Dialog.addMessage("No 'aaa000BEADS.tif' image found");
	Dialog.addMessage("Code can't work without bead image stack saved in this folder");
	Dialog.show();
}

selectWindow(beadImgMM); nsBead=nSlices; close();

if (nf!=(nsBead+1))
{
	Dialog.create("WARNING");
	Dialog.addMessage("Mismatch between number of bead images and number of images to be processed for anisotropy");
	Dialog.addMessage("If 2 images can be registered with 1 bead image, duplicate the bead image sequentially in the stack");
	Dialog.show();
}


for (i=1; i<=nsBead; i++)
{
open(dir+list[beadImgIndex]);
selectWindow(beadImgMM); setSlice(i); run("Duplicate...", "title=beads");
selectWindow(beadImgMM); close();

beadIMAGEfilename="beads";
selectWindow(beadIMAGEfilename); cellWid=getWidth(); cellHgt=getHeight(); //images of beads and cells should be of same size

open(dir+list[i]);
openimg=nImages;
for (imgID=0; imgID<openimg; imgID++)
{
	selectImage(imgID+1);
	cellIMAGEfilename=getTitle();
	
	if ((cellIMAGEfilename!=beadIMAGEfilename)&&(cellIMAGEfilename!=GfactorImage))
	{
		selectWindow(cellIMAGEfilename);
		slices=nSlices;//number of slices in image stack of cells
//New Images to store Anisotropy and Total Intensity of the whole field of view
		newImage("WholeFieldAni", "32-bit black", round(cellWid/2), cellHgt, slices);
		newImage("WholeFieldInt", "32-bit black", round(cellWid/2), cellHgt, slices);
		selectWindow(cellIMAGEfilename);
		
//divide image into small boxes (to perform individual registration on them, for best results), and add selections to ROI Manager
		run("ROI Manager...");
		nROI=roiManager("count");
		if (nROI>0) { roiManager("Deselect"); roiManager("Delete"); }
//we want to divide a 1024X2048 image into 5 overlapping rows and 3 overlapping columns
		NboxRow=4; NboxCol=2; rboxsz=round(cellHgt/NboxRow); xxr=0; yyr=0;
		for (iro=0; iro<(NboxRow+1); iro++) //to fit 512 five times in 2048, we overlap 512/4=128 sized sections
		{
			for (ico=0; ico<(NboxCol+1); ico++) //to fit 512 three times in 1024, we overlap 512/2=256 sized sections
			{
				xxr=ico*rboxsz/NboxCol; yyr=iro*(rboxsz-(rboxsz/NboxRow));
				makeRectangle(xxr, yyr, rboxsz, rboxsz); roiManager("Add");
			}
		}
//select an ROI, perform Anisotropy calculation
nROI=roiManager("count"); //find number of ROIs stored
//Remove slice number from image
for (iROI=0; iROI<nROI; iROI++)
{
	roiManager("Select", iROI);
	roiManager("Remove Slice Info");
}

for (iROI=0; iROI<nROI; iROI++)
{
	
//create rectangular ROI from polygon selection for registration
	selectWindow(cellIMAGEfilename);
	roiManager("Select", iROI);
	getSelectionCoordinates(Xp, Yp);
	x1=Xp[0]; y1=Yp[0]; w=abs(Xp[0]-Xp[1]); h=abs(Yp[1]-Yp[2]);
	x2=x1+xshift; y2=y1+yshift;

//add suffix after any number in cell image name, or at the end if numbers are absent
	lc=lengthOf(cellIMAGEfilename);
//remove any file name extension
	if (substring(cellIMAGEfilename,(lc-4),(lc-3))==".") { cellIMAGEname=substring(cellIMAGEfilename,0,(lc-4)); lc=lc-4; }
	else { cellIMAGEname=cellIMAGEfilename; }
	
//use Descriptor-based Registration Plugin to register the image subset of beads in the ROI
	selectWindow(beadIMAGEfilename);
	makeRectangle(x1, y1, w, h);
	run("Duplicate...", "title=Bperpendicular");
	selectWindow("Bperpendicular");
	selectWindow(beadIMAGEfilename);
	makeRectangle(x2, y2, w, h);
	run("Duplicate...", "title=Bparallel");
	selectWindow("Bparallel");
	run("Multiply...", "value="+gFactorUNPolarised);
	selectWindow("Bparallel");
	//run("Merge Channels...", "c1=Bparallel c2=Bperpendicular create keep ignore");
	//run("Descriptor-based registration (2d/3d)", "first_image=Bperpendicular second_image=Bparallel brightness_of=[Interactive ...] approximate_size=[Interactive ...] type_of_detections=[Interactive ...] subpixel_localization=[3-dimensional quadratic fit] transformation_model=["+TransformModel+" (2d)] images_pre-alignemnt=[Not prealigned] number_of_neighbors=3 redundancy=1 significance=3 allowed_error_for_ransac=5 choose_registration_channel_for_image_1=1 choose_registration_channel_for_image_2=1 create_overlayed add_point_rois");
	run("Descriptor-based registration (2d/3d)", "first_image=Bperpendicular second_image=Bparallel brightness_of=Medium approximate_size=["+beadsz+" px] type_of_detections=[Maxima only] subpixel_localization=[3-dimensional quadratic fit] transformation_model=["+TransformModel+" (2d)] images_pre-alignemnt=[Not prealigned] number_of_neighbors=3 redundancy=1 significance=3 allowed_error_for_ransac=5 choose_registration_channel_for_image_1=1 choose_registration_channel_for_image_2=1 create_overlayed add_point_rois");
	selectWindow("Bparallel");	close();	
	selectWindow("Bperpendicular");	close();
	selectWindow("Fused Bperpendicular & Bparallel"); close();

//now use Descriptor-based Registration Plugin to register the Gfactor image subset using the same transformation matrix as used for beads
	if (gfm==1)
	{
		selectWindow(GfactorImage);
		makeRectangle(x1, y1, w, h);
		run("Duplicate...", "title=GFperpendicular");
		selectWindow("GFperpendicular");
		selectWindow(GfactorImage);
		makeRectangle(x2, y2, w, h);
		run("Duplicate...", "title=GFparallel");
		selectWindow("GFparallel");
		run("Descriptor-based registration (2d/3d)", "first_image=GFperpendicular second_image=GFparallel reapply");
		selectWindow("GFparallel");	close();
		selectWindow("GFperpendicular"); close();
		selectWindow("Fused GFperpendicular & GFparallel");
		run("Hyperstack to Stack");
		selectWindow("Fused GFperpendicular & GFparallel");
		run("Stack to Images");
		imageCalculator("Divide create 32-bit", "GFparallel","GFperpendicular");
		selectWindow("GFparallel");	close();
		selectWindow("GFperpendicular"); close();
		selectWindow("Result of GFparallel"); rename("Gf");
	}

//now use Descriptor-based Registration Plugin to register the cell image subset using the same transformation matrix as used for beads
	selectWindow(cellIMAGEfilename);
	makeRectangle(x1, y1, w, h);
	run("Duplicate...", "title=Cperpendicular duplicate");	selectWindow("Cperpendicular");	
	selectWindow(cellIMAGEfilename);
	makeRectangle(x2, y2, w, h);
	run("Duplicate...", "title=Cparallel duplicate");	selectWindow("Cparallel");
	//run("Merge Channels...", "c1=Cparallel c2=Cperpendicular create keep ignore");
	run("Descriptor-based registration (2d/3d)", "first_image=Cperpendicular second_image=Cparallel reapply");
	selectWindow("Fused Cperpendicular & Cparallel");
	selectWindow("Cparallel"); close();
	selectWindow("Cperpendicular"); close();
//put Parallel Channel as Channel 1 and Perpendicular Channel as Channel 2 in a new Hyperstack
	selectWindow("Fused Cperpendicular & Cparallel");
	run("Duplicate...", "title=Cperpendicular duplicate channels=1");
	selectWindow("Fused Cperpendicular & Cparallel");
	run("Duplicate...", "title=Cparallel duplicate channels=2");
	selectWindow("Fused Cperpendicular & Cparallel"); close();
	run("Concatenate...", "  title=[Fused Cparallel & Cperpendicular] open image1=Cparallel image2=Cperpendicular image3=[-- None --]");
	selectWindow("Fused Cparallel & Cperpendicular");
	run("Stack to Hyperstack...", "order=xytcz channels=2 slices=1 frames="+slices+" display=Color");
	regIMAGEfilename="Fused Cparallel & Cperpendicular";
	
	if (pflag==1)
	{
		selectWindow(regIMAGEfilename);
		run("Duplicate...", "title=registeredImage duplicate");
		selectWindow("registeredImage");
		regIMAGEfilenameSave=cellIMAGEname+"-registered.tif";
		saveAs("Tiff", savedfilepath+regIMAGEfilenameSave); close();
	}

//end of image-registration module
//PERFORM ANISOTROPY CALCULATION

	selectWindow(regIMAGEfilename);	setSlice(1);

	if (StRegflag==1) //to further cancel any registration disputes (backup)
	{
		run("Linear Stack Alignment with SIFT", "initial_gaussian_blur=1.60 steps_per_scale_octave=3 minimum_image_size=64 maximum_image_size=2048 feature_descriptor_size=4 feature_descriptor_orientation_bins=8 closest/next_closest_ratio=0.92 maximal_alignment_error=25 inlier_ratio=0.05 expected_transformation=Affine interpolate");
		selectWindow(regIMAGEfilename); close();
		selectWindow("Aligned "+(slices*2)+" of "+(slices*2));
		rename(regIMAGEfilename);
	}

	selectWindow(regIMAGEfilename);
	run("Duplicate...", "title=numP duplicate channels=1"); selectWindow("numP"); run("32-bit"); //parallel channel
	selectWindow(regIMAGEfilename);
	run("Duplicate...", "title=numS duplicate channels=2"); selectWindow("numS"); run("32-bit"); //perpendicular channel
	selectWindow(regIMAGEfilename); close();

//Subtract respective background values
	selectWindow("numP"); run("Subtract...", "value="+backParallel+" stack");
	selectWindow("numS"); run("Subtract...", "value="+backPerpendicular+" stack");
	
//Multiply G-factor with perpendicular channel
	if (gfm==1) { imageCalculator("Multiply 32-bit stack", "numS","Gf"); selectWindow("Gf"); close(); }
	if (gfm==0) { run("Multiply...", "value="+Gfactor+" stack"); }
	
	selectWindow("numP");	run("Duplicate...", "title=denP duplicate");
	selectWindow("numS");	run("Duplicate...", "title=denS duplicate");
	selectWindow("denS");   run("Multiply...", "value=2 stack");
	imageCalculator("Subtract create 32-bit stack", "numP","numS"); selectWindow("numP"); close(); selectWindow("numS"); close();
	selectWindow("Result of numP");
	imageCalculator("Add create 32-bit stack", "denP","denS"); selectWindow("denP"); close(); selectWindow("denS"); close();
	selectWindow("Result of denP"); rename("TotalIntensity"); //Total-Intensity Image
	imageCalculator("Divide create 32-bit stack", "Result of numP","TotalIntensity"); selectWindow("Result of numP"); close();
	selectWindow("Result of Result of numP");  rename("Anisotropy"); //Anisotropy Image

//create mask of cell image based on intensity values
	if (maskcheck==1)
	{
		selectWindow("TotalIntensity"); run("Grays");
		run("Duplicate...", "title=mask duplicate");
		//selectWindow("mask"); run("Duplicate...", "title=maskCheck duplicate");
		selectWindow("mask");
		//setAutoThreshold("Huang dark");
		//setAutoThreshold("Triangle dark");
		setAutoThreshold("Percentile dark");
		///*
		run("NaN Background", "stack");
		run("Multiply...", "value=0 stack");
		run("Add...", "value=1 stack");
		//*/
		/*
		run("Convert to Mask", "method=Triangle background=Dark calculate black");
		selectWindow("mask"); run("Divide...", "value=255.000 stack"); run("32-bit");
		*/
		imageCalculator("Multiply 32-bit stack","TotalIntensity","mask");
		imageCalculator("Multiply 32-bit stack","Anisotropy","mask");
		selectWindow("mask"); close();
	}

//SEQUENTIALLY copy Anisotropy and Total Intensity image slices and paste them into WholeFieldAni and WholeFieldInt per slice per ROI
	for (iimage=1; iimage<=slices; iimage++)
	{
		//copy Anisotropy image slices into WholeFieldAni as per ROI
		selectWindow("Anisotropy"); setSlice(iimage);
		makeRectangle(xCutOff, yCutOff, w-(2*xCutOff), h-(2*yCutOff)); run("Copy");
		selectWindow("WholeFieldAni"); setSlice(iimage);
		roiManager("Select", iROI);
		getSelectionCoordinates(Xpp, Ypp); nppp=Xpp.length; Array.getStatistics(Xpp, xppmin, xppmax, xppmean, xppstd); Array.getStatistics(Ypp, yppmin, yppmax, yppmean, yppstd);
		makeRectangle(xppmin+xCutOff, yppmin+yCutOff, w-(2*xCutOff), h-(2*yCutOff));
		selectWindow("WholeFieldAni"); run("Paste");
		//copy Total Intensity image slices into WholeFieldAni as per ROI
		selectWindow("TotalIntensity"); setSlice(iimage);
		makeRectangle(xCutOff, yCutOff, w-(2*xCutOff), h-(2*yCutOff)); run("Copy");
		selectWindow("WholeFieldInt"); setSlice(iimage);
		roiManager("Select", iROI);
		getSelectionCoordinates(Xpp, Ypp); nppp=Xpp.length; Array.getStatistics(Xpp, xppmin, xppmax, xppmean, xppstd); Array.getStatistics(Ypp, yppmin, yppmax, yppmean, yppstd);
		makeRectangle(xppmin+xCutOff, yppmin+yCutOff, w-(2*xCutOff), h-(2*yCutOff));
		selectWindow("WholeFieldInt"); run("Paste");
	}
	selectWindow("Anisotropy"); close();
	selectWindow("TotalIntensity"); close();
	
}// end of- for (iROI=0; iROI<nROI; iROI++)

//SAVE RESULTS AFTER STITCHING BACK THE BOXES
	selectWindow("WholeFieldInt");	run("Select None"); run("Enhance Contrast", "saturated=0.35"); run("Grays");
	saveAs("Tiff", savedfilepath+"TotalIntensity32bit-"+cellIMAGEname); close();
	selectWindow("WholeFieldAni"); run("Select None"); setMinAndMax(minAni, maxAni); run("16 colors");
	saveAs("Tiff", savedfilepath+"Anisotropy32bit-"+cellIMAGEname); close();

	
	}// end of- if ((cellIMAGEfilename!=beadIMAGEfilename)&&(cellIMAGEfilename!=GfactorImage))
}// end of- for (imgID=0; imgID<openimg; imgID++)

roiManager("Deselect"); roiManager("Delete");
selectWindow(cellIMAGEfilename); close();
selectWindow(beadIMAGEfilename); close();
}