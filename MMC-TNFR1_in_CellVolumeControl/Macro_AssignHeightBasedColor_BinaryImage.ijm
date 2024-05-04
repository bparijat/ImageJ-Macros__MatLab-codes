//Macro to assign a specific height-dependent colour to a slice in a Z stack BINARY IMAGE
//and create a final pseudo coloured image that indicates all the heights in specific colour
//best to have the voxel size in-built into the image

maxHeight=12; //in microns, for normalizing colors
makeColorBar=0; //0=false,1=true
reverseLUT=1; //0=false,1=true

MinPixelValue=0; //Minimum grayvalue in image to be put in colorbar scale
MaxPixelValue=maxHeight; //Maximum grayvalue in image to be put in colorbar scale
decimalPlaces=2;
LUTtoApply="royal"; //"phase";

id=getTitle; getVoxelSize(width, height, depth, unit); ns=nSlices; //get details of image stack
finalName=id+"_ColoredByHeight";
selectImage(id); run("Select None"); run("Duplicate...", "title=C duplicate");
selectImage("C"); run("32-bit"); run("Divide...", "value=255 stack"); //convert to 32 bit image and make all 255's 1 for multiplication with height values
selectImage("C"); run("Select None"); run("Duplicate...", "title=H duplicate"); //for 2D height map

//Multiply each slice with their height value above the basal plane for assigning height-baed color
for (i=1; i<=ns; i++)
{
	selectImage("C"); setSlice(i); run("Multiply...", "value="+(depth*(i-1))+" slice");
}
selectImage("C"); setMinAndMax(MinPixelValue, MaxPixelValue);
//Assign colors based on height values
selectImage("C"); run(LUTtoApply); getLut(R, G, B);
//Reverse LUT for better visibility
if (reverseLUT==1)
{
	Array.reverse(R);
	Array.reverse(G);
	Array.reverse(B);
	R[0]=0; B[0]=0; G[0]=0;
	selectImage("C"); setLut(R, G, B);
}

selectImage("C"); run("Select None"); run("Duplicate...", "title="+finalName+" duplicate"); setSlice(1); close("C");

//For 2D height map, multiply each slice with their depth value
for (i=1; i<=ns; i++)
{
	selectImage("H"); setSlice(i); run("Multiply...", "value="+(depth)+" slice");
}
selectImage("H"); run("Z Project...", "projection=[Sum Slices]"); close("H");
selectImage("SUM_H");  setMinAndMax(MinPixelValue, MaxPixelValue); setLut(R, G, B);
selectImage("SUM_H"); run("Select None"); rename(finalName+"_HeightMap2D");

//Create color bar depicting exact height
if (makeColorBar==1)
{
	newImage("CalibrationBar", "32-bit white", 2048, 2048, 1);
	run("Add...", "value=99999999999");
	setMinAndMax(MinPixelValue, MaxPixelValue); setLut(R, G, B);
	run("Calibration Bar...", "location=[Separate Image] fill=[White] label=Black number=5 decimal="+decimalPlaces+" font=10 zoom=15");
	close("CalibrationBar"); selectWindow("CBar"); rename("CBar_HeightMicrons");
}