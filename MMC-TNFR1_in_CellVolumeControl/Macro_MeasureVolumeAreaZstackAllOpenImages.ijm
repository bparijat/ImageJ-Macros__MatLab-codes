//This Macro measures the volume of an object in a Z stack, for ALL OPEN IMAGES
//Auto-threshold with NaN background, convert foreground to 1, count number of foreground pixels, multiply with pixel size parameters
//CLICK ON IMAGE TO GET TITLE OF IMAGE then RUN macro

nOpenImgs=nImages;
for (imgID1=0; imgID1<nOpenImgs; imgID1++)
{
selectImage(imgID1+1); fname=getTitle(); run("32-bit"); slices=nSlices;
getVoxelSize(width, height, depth, unit);

//setAutoThreshold("Default dark stack");
setAutoThreshold("Huang dark stack");
//setAutoThreshold("Otsu dark stack");
run("NaN Background", "stack");
run("Multiply...", "value=0 stack");
run("Add...", "value=1 stack");
run("Set Measurements...", "integrated redirect=None decimal=3");
run("Clear Results");
IntValA=0;
selectWindow(fname); setSlice(1); run("Measure");
ka=parseInt(getResultString("RawIntDen", (0))); //basal area measurement
if (ka>0) { IntValA=IntValA+(ka*width*height); }
run("Clear Results");
IntVal=0;
for (i=0; i<slices; i++)
{
	selectWindow(fname); setSlice(i+1); run("Measure");
	k=parseInt(getResultString("RawIntDen", (i)));
	if (k>0) { IntVal=IntVal+(k*width*height*depth); }
}
print(fname+"\t"+"volume"+"="+"\t"+IntVal+"\t"+"area"+"="+"\t"+IntValA);
selectWindow(fname); //run("Multiply...", "value=255 stack"); setSlice(1); setMinAndMax(0, 255);
}
run("Close All");