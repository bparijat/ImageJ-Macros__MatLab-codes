/*
 * opens all files in the directory and calculates intensity and anisotropy
 * within intensity or anisotropy based masks made from 32-bit total intensity images
 */

IntStartAlpha="T"; AniStartAlpha="A";
threshmin=9000;  threshmax=140000;
minAniTh=0.11; maxAniTh=0.35;
Dialog.create("Operation Details");
Dialog.addCheckbox("Use Auto-Threshold", true);
Dialog.addMessage("For User-defined Thresholding, uncheck Auto-Threshold box and input values here");
Dialog.addNumber("Minimum Intensity Threshold:", threshmin);
Dialog.addNumber("Maximum Intensity Threshold:", threshmax);
Dialog.addMessage("Check box to additionally segment using anisotropy values");
Dialog.addCheckbox("Use Anisotropy Values", false);
Dialog.addNumber("Minimum Anisotrpy Threshold:", minAniTh);
Dialog.addNumber("Maximum Anisotrpy Threshold:", maxAniTh);
Dialog.show();
chkbx0=Dialog.getCheckbox();
chkbx1=Dialog.getCheckbox();
dir = getDirectory("Browse to Folder");
list = getFileList(dir);
nf=list.length;

//FIND start and end number of Anisotropy and Intensity files
stA=1; stT=1; enA=1; enT=1;
for (i=0; i<nf; i++)  { if (startsWith(list[i], AniStartAlpha)) { stA=i; break; } }
for (i=0; i<nf; i++)  { if (startsWith(list[i], IntStartAlpha)) { stT=i; break; } }
for (i=0; i<nf; i++)
{
	if (startsWith(list[i], AniStartAlpha)) { enA=i; }
	if (startsWith(list[i], IntStartAlpha)) { enT=i; }
}
nCells=enA-stA+1; //number of Anisotropy images

//OBTAIN typical image details in the containing folder
ns00=0; maxslice=0;
for (i=stA; i<=enA; i++)
{
	open(dir+list[i+1]);
	AniImg=getTitle();
	ns00=nSlices;//number of slices in image stack
	selectWindow(AniImg);
	close();
	if (maxslice<ns00)  { maxslice=ns00; }
}

//CREATE an Array to store measured values in Result table
IntVal=newArray(maxslice*nCells);  AniVal=newArray(maxslice*nCells);  AreaVal=newArray(maxslice*nCells);

//OPEN Anisotropy and Intensity images
//MAKE masks from intensity images and measure the anisotropy values within it

run("Set Measurements...", "area mean modal min center median redirect=None decimal=3");

for (iA=stA,iT=stT; iA<=enA; iA++,iT++)
{
	open(dir+list[iA]);
	AniImg=getTitle();
	run("Select All");
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	open(dir+list[iT]);
	IntImg=getTitle();
	run("Select All");
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	
	selectWindow(AniImg);
	setMinAndMax(0.0800, 0.2500);
	selectWindow(IntImg);
	run("Enhance Contrast", "saturated=0.35");
	nslice1=nSlices;//number of slices in current image stack
	
	//CREATE masks from Intensity image slices and add them to ROI manager
	run("ROI Manager...");
	nroi=roiManager("count");
	if (nroi>0) { roiManager("Deselect"); roiManager("Delete"); }
	//make a mask based on intensity
	selectWindow(IntImg);  setSlice(1);  run("Select All");
	run("Duplicate...", "title=Imask duplicate");
	selectWindow("Imask");  nslice2=nSlices;  Imw=getWidth();  Imh=getHeight();
	if (chkbx0==true)	{ for (islice2=1; islice2<=nslice2; islice2++)	 { setSlice(islice2); setAutoThreshold("Default dark"); run("NaN Background", "slice"); }	}
	if (chkbx0==false) { setThreshold(threshmin, threshmax); run("NaN Background", "stack"); }
	run("Multiply...", "value=0.000 stack");
	run("Add...", "value=1.000 stack");
	newImage("Imask2", "32-bit black", Imw, Imh, nslice2);
	run("Clear Results");
	selectWindow("Imask");
	for (islice2=0; islice2<nslice2; islice2++)
	{
		selectWindow("Imask"); setSlice(islice2+1); run("Measure");
		cX=parseFloat(getResultString("XM", islice2)); cY=parseFloat(getResultString("YM", islice2));
		doWand(cX, cY); roiManager("Add");
		selectWindow("Imask2"); roiManager("Select", islice2);
		run("Add...", "value=1 slice");
		setAutoThreshold("Default dark"); run("NaN Background", "slice");
	}
	selectWindow("Imask"); close();
	selectWindow("Imask2"); rename("Imask");
	//multiply Intensity and Anisotropy images with Imask
	selectWindow(IntImg);  setSlice(1);
	imageCalculator("Multiply 32-bit stack", IntImg, "Imask");
	selectWindow(AniImg);  setSlice(1);
	imageCalculator("Multiply 32-bit stack", AniImg, "Imask");
	selectWindow("Imask"); close();
	//segmentation based on anisotropy values
	if (chkbx1==true)
	{
		selectWindow(AniImg);  setSlice(1);  run("Select All");
		//SET threshold intensities for 32 bit image
		run("Duplicate...", "title=Amask duplicate");
		selectWindow("Amask");
		setThreshold(minAniTh, maxAniTh);
		run("NaN Background", "stack");
		run("Multiply...", "value=0.000 stack");
		run("Add...", "value=1.000 stack");
		//multiply Intensity and Anisotropy images with Imask
		selectWindow(IntImg);  setSlice(1);
		imageCalculator("Multiply 32-bit stack", IntImg, "Amask");
		selectWindow(AniImg);  setSlice(1);
		imageCalculator("Multiply 32-bit stack", AniImg, "Amask");
		selectWindow("Amask"); //close();
	}
	//MEASURE Intensity and Anisotropy values in Thresholded images
	run("Clear Results");
	for (j=0; j<nslice1; j++) {	selectWindow(IntImg); roiManager("Select", j); run("Measure"); mesI=parseFloat(getResultString("Mean", j)); mesA=parseFloat(getResultString("Area", j)); IntVal[(iA*maxslice+j)]=mesI; AreaVal[(iA*maxslice+j)]=mesA; }
	//MEASURE and STORE Anisotropy values
	run("Clear Results");
	for (j=0; j<nslice1; j++) { selectWindow(AniImg); roiManager("Select", j); run("Measure"); mesA=parseFloat(getResultString("Mean", j)); AniVal[(iA*maxslice+j)]=mesA; }
	
	selectWindow(AniImg);	close();
	selectWindow(IntImg);	close();
	roiManager("Deselect"); roiManager("Delete");
}

//PRINT measured values
///*
strLine=" ";
for (i=0; i<nCells; i++)
{
	strLine=strLine+"Int_"+(i+1)+"\t";
}
print(strLine);
strLine=" ";
for (j=0; j<maxslice; j++)
{
	for (i=0; i<nCells; i++)
	{
		strLine=strLine+IntVal[(i*maxslice+j)]+"\t";
	}
	print(strLine);
	strLine=" ";
}
print("\n \n");
strLine=" ";
for (i=0; i<nCells; i++)
{
	strLine=strLine+"Ani_"+(i+1)+"\t";
}
print(strLine);
strLine=" ";
for (j=0; j<maxslice; j++)
{
	for (i=0; i<nCells; i++)
	{
		strLine=strLine+AniVal[(i*maxslice+j)]+"\t";
	}
	print(strLine);
	strLine=" ";
}
print("\n \n");
strLine=" ";
for (i=0; i<nCells; i++)
{
	strLine=strLine+"Area_"+(i+1)+"\t";
}
print(strLine);
strLine=" ";
for (j=0; j<maxslice; j++)
{
	for (i=0; i<nCells; i++)
	{
		strLine=strLine+AreaVal[(i*maxslice+j)]+"\t";
	}
	print(strLine);
	strLine=" ";
}
//*/

//PRINT measured values in the format to create box-plot
/*
print("Intensity");
strLine=" ";
for (i=0; i<nCells; i++)
{
	for (j=0; j<maxslice; j++)
	{
		strLine=strLine+IntVal[(i*maxslice+j)]+"\t";
	}
	print(strLine);
	strLine=" ";
}
print("\n \n");
print("Anisotropy");
strLine=" ";
for (i=0; i<nCells; i++)
{
	for (j=0; j<maxslice; j++)
	{
		strLine=strLine+AniVal[(i*maxslice+j)]+"\t";
	}
	print(strLine);
	strLine=" ";
}
*/