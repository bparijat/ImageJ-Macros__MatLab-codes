fname=getTitle();
nROI=roiManager("count");
for (i=0; i<nROI; i++)
{
	selectWindow(fname);
	roiManager("Select", i);
	run("Duplicate...", "title=["+fname+"-"+i+1+"] duplicate");
	run("Make Inverse");
	run("Multiply...", "value=0 stack");
	run("Select None");
}