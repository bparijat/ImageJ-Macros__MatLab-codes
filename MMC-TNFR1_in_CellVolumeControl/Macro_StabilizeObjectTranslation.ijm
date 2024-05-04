//MACRO to Stabilize Translation of an Object in a video
//we will align all slices in a stack according to center of mass of an object whose ROIs are known
//click on IMAGE STACK you want to align, to get its title, then run macro

imname=getTitle(); ns=nSlices; w=getWidth; h=getHeight; hw=round(w/2); bd=bitDepth();

nROI=roiManager("count"); //find number of ROIs stored

CMx=newArray(nROI); CMy=newArray(nROI);
Xshift=newArray(nROI); Yshift=newArray(nROI);
//find max shift between positions
maxShiftX=0; maxShiftY=0;
for (i = 0; i < nROI; i++)
{
	roiManager("Select", i);
	List.setMeasurements; //(List.setMeasurements can be used for any selection area, doesn't have to be in ROI Manager)
	CMx[i]=List.getValue("XM"); CMy[i]=List.getValue("YM"); //get Object's Center of Mass in slice i
	Xshift[i]=round(CMx[i]-CMx[0]); Yshift[i]=round(CMy[i]-CMy[0]); //find shift in pixels
	if (abs(Xshift[i])>maxShiftX) { maxShiftX=abs(Xshift[i]); }
	if (abs(Yshift[i])>maxShiftY) { maxShiftY=abs(Yshift[i]); }
}

//make new image stack to adjust positions
newImage("egg", bd+"-bit black", w+2*maxShiftX, h+2*maxShiftY, ns);

//to shift slice i and align with the first slice
for (i = 0; i < nROI; i++)
{
	selectWindow(imname); roiManager("Select", i); run("Select All"); run("Copy");
	selectWindow("egg"); setSlice(i+1); makeRectangle(maxShiftX-Xshift[i], maxShiftY-Yshift[i], w, h); run("Paste");
}
selectWindow(imname); run("Select None"); setSlice(1);
selectWindow("egg"); rename(imname+"-Translation_Stabilised"); setSlice(1); run("Select None"); run("Enhance Contrast", "saturated=0.35");