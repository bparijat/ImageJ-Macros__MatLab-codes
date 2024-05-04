//MACRO TO CALCULATE NUCLEUS/CYTOPLASMIC RATIO USING RAW INTENSITY VALUES
//save all ROIs for nuclei and cell body in different zip files
//open the image/stack for analysis, click on the image, run macro and follow prompts
//the results will be copied to clipboard, paste them in a spreadsheet

imname=getTitle(); //get name of image to be analyzed
roiManager("deselect"); roiManager("delete");

//find path of ROIset files
pathNucROI = File.openDialog("Select ROI set for nuclei");
pathCytROI = File.openDialog("Select ROI set for cell");

//collect ROI info of all the nuclei
roiManager("open", pathNucROI);
nROInuc=roiManager("count");
//store centroids of all ROIs
Xnuc=newArray(nROInuc); Ynuc=newArray(nROInuc); AreaNuc=newArray(nROInuc);
for (i=0; i<=(nROInuc-1); i++)
{
	selectWindow(imname); roiManager("Select", i);
	Xnuc[i]=getValue("X"); Ynuc[i]=getValue("Y");  AreaNuc[i]=getValue("Area");
}
roiManager("deselect"); roiManager("delete");

//collect ROI info of all the cell bodies
roiManager("open", pathCytROI);
nROIcyt=roiManager("count");
//store centroids of all ROIs
Xcyt=newArray(nROIcyt); Ycyt=newArray(nROIcyt); AreaCyt=newArray(nROIcyt);
for (i=0; i<=(nROIcyt-1); i++)
{
	selectWindow(imname); roiManager("Select", i);
	Xcyt[i]=getValue("X"); Ycyt[i]=getValue("Y"); AreaCyt[i]=getValue("Area");
}
roiManager("deselect"); roiManager("delete");

//identify the nucleus ROI for the specific cell ROI by finding minimum distance between their centroids
Array.getStatistics(AreaCyt, minCellArea, maxCellArea, meanCellArea, stdDevCellArea);
CytNucMATCH=newArray(nROIcyt); BadCellROI=newArray(nROIcyt);
for (i=0; i<=(nROIcyt-1); i++)
{
	minDist=99999999999;
	BadCellROI[i]=0; //if no corresponding nucleus ROI found, ROI is marked as a bad cell (=1)
	for (j=0; j<=(nROInuc-1); j++)
	{
		xDIFF=pow((Xcyt[i]-Xnuc[j]),2); yDIFF=pow((Ycyt[i]-Ynuc[j]),2); dist=sqrt(xDIFF+yDIFF);
		if (dist<=minDist) { minDist=dist; CytNucMATCH[i]=j; } //identify corresponding nucleus ROI
	}
	if (minDist>sqrt(maxCellArea)) { BadCellROI[i]=1; } //mark bad cells
	print(i,"\t",CytNucMATCH[i],"\t",minDist,"\t",BadCellROI[i],"\t",Xnuc[i],"\t",Xcyt[i],"\t",Ynuc[i],"\t",Ycyt[i]);
}

//open cell ROIs, measure intensities in cell and corresponding nucleus ROIs, calculate ratios
IntensityCELL=newArray(nROIcyt); IntensityNUC=newArray(nROIcyt); IntensityCYT=newArray(nROIcyt);
RatioNucCell=newArray(nROIcyt); RatioNucCyt=newArray(nROIcyt);
print("Cell Area","\t","Nucleus Area","\t","Whole Cell Intensity","\t","Nucleus Intensity","\t","Cytoplasm Intensity","\t","Nucleus/WholeCell","\t","Nucleus/Cytoplasm");
for (i=0; i<=(nROIcyt-1); i++)
{
	if (BadCellROI[i]==0) //check if cell is good
	{
		//measure whole cell's total intensity
		roiManager("open", pathCytROI);
		selectWindow(imname); roiManager("Select", i);
		IntensityCELL[i]=getValue("RawIntDen");
		roiManager("deselect"); roiManager("delete");
		//measure nucleus' total intensity
		roiManager("open", pathNucROI);
		selectWindow(imname); roiManager("Select", CytNucMATCH[i]);
		IntensityNUC[i]=getValue("RawIntDen");
		roiManager("deselect"); roiManager("delete");
		//calculate ratios
		IntensityCYT[i]=abs(IntensityCELL[i]-IntensityNUC[i]);
		RatioNucCell[i]=IntensityNUC[i]/IntensityCELL[i];
		RatioNucCyt[i]=IntensityNUC[i]/IntensityCYT[i];
		//print(AreaCyt[i],"\t",AreaNuc[i],"\t",IntensityCELL[i],"\t",IntensityNUC[i],"\t",IntensityCYT[i],"\t",RatioNucCell[i],"\t",RatioNucCyt[i]);
	}
}
/*
selectWindow(imname); close();
AllData=getInfo("log"); String.copy(AllData); //copy everything in Log window
Dialog.create("MESSAGE");
Dialog.addMessage("All Results copied");
Dialog.addMessage("just PASTE them in a spreadsheet, then click on OK");
Dialog.addMessage("(clicking on CANCEL stops the program)");
Dialog.setLocation(0,0)
Dialog.show();
print("\\Clear"); //clear Log window
*/