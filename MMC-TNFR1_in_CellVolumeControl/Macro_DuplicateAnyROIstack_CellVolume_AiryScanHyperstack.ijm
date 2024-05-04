fname=getTitle();
getDimensions(width, height, channels, slices, frames);

nROI=roiManager("count");
for (i=0; i<nROI; i++)
{
	selectWindow(fname);
	roiManager("Select", i);
	Roi.getPosition(channelN, sliceN, frameN);
	run("Duplicate...", "title=["+fname+"-"+i+1+"] duplicate frames="+frameN);
	run("Make Inverse");
	run("Multiply...", "value=0 stack");
	run("Add...", "value=10003 stack");
	run("Add Specified Noise...", "stack standard=1.28");
	run("Select None");
}