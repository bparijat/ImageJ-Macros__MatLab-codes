//Macro to calculate time resolved anisotropy between the parallel channel 'P' and the perpendicular channel 'S'
//This macro includes additional smoothing between frames
//Anisotropy r = (I_p - Gfactor*I_s)/(I_p + 2*Gfactor*I_s)
//create duplicates of parallel channel named 'P' and perpendicular channel named 'S'
//run macro

Gfactor=3.012; //from fluorescein-in-water measurements (multiply with s-channel)
background=0; //background count
intTime=50; //in ns, represents total integration time
binsz=3; //in pixels, indicating radius of pixel binning
minAni=0.000000000;
maxAni=0.400000000;
smoothT=1; //number of time frames to be used for applying moving average

//GET USER INPUTS
Dialog.create("Operation Details");
Dialog.addMessage("SEE THAT THE PARALLEL CHANNEL HAS BEEN RENAMED TO 'P'\nAND THE PERPENDICULAR CHANNEL HAS BEEN RENAMED TO 'S'\n");
Dialog.addNumber("Use this G-factor value:", Gfactor);
Dialog.addNumber("Average Background (parallel channel)", background);
Dialog.addNumber("Total integration time (in ns):", intTime);
Dialog.addNumber("Minimum Anisotrpy displayed:", minAni);
Dialog.addNumber("Maximum Anisotrpy displayed:", maxAni);
Dialog.addMessage("To smooth images (over X,Y) before Anisotropy calculations by binning pixels...");
Dialog.addChoice("...Pixel binning method:", newArray("None","Sum","Average"));
Dialog.addNumber("use n-by-n pixel binning (n MUST BE AN ODD NUMBER), where n=", binsz);
Dialog.addMessage("To smooth pixels over TIME before Anisotropy calculations...");
Dialog.addMessage("...by computing moving average over time...");
Dialog.addNumber("...use T_avg=", smoothT);
Dialog.addMessage("Please use an ODD NUMBER for T_avg  (To not apply moving average, use T_avg=1)");
Dialog.addCheckbox("Tick checkbox to keep smoothed original images after running macro", false);
Dialog.show();
Gfactor=Dialog.getNumber();
background=Dialog.getNumber();
intTime=Dialog.getNumber();
minAni=Dialog.getNumber();
maxAni=Dialog.getNumber();
BinModel=Dialog.getChoice();
binsz=Dialog.getNumber();
smoothT=Dialog.getNumber();
keepSmooth=Dialog.getCheckbox();

while ((binsz%2)==0) { Dialog.create("ERROR!"); Dialog.addNumber("Binning cannot be performed with n = even number\n, input n = odd number", binsz); Dialog.show(); binsz=Dialog.getNumber(); }
while ((smoothT%2)==0) { Dialog.create("ERROR!"); Dialog.addNumber("Binning cannot be performed with T_avg = even number\n, input T_avg = odd number", smoothT); Dialog.show(); smoothT=Dialog.getNumber(); }

//GENERATE BINNING KERNELS ALONG X and Y DIRECTIONS
kernelBinXY="";
if (BinModel=="None") {	kernelBinXY="text1=[0 0 0\n0 1 0\n0 0 0] stack"; }
if (BinModel=="Sum")
{
	for (i=1; i<=binsz; i++)
	{
		for (j=1; j<=binsz; j++)
		{
			kernelBinXY=kernelBinXY+"1 ";
		}
		kernelBinXY=kernelBinXY+"\n";
	}
	kernelBinXY="text1=["+kernelBinXY+"] stack";
}
if (BinModel=="Average")
{
	for (i=1; i<=binsz; i++)
	{
		for (j=1; j<=binsz; j++)
		{
			kernelBinXY=kernelBinXY+"1 ";
		}
		kernelBinXY=kernelBinXY+"\n";
	}
	kernelBinXY="text1=["+kernelBinXY+"] normalize stack";
}

//GENERATE BINNING KERNELS ALONG T(Z) DIRECTION FOR APPLYING MOVING AVERAGE
kernelBinT="";
spinalPos=floor(smoothT/2)+1;
centralPixel=1/smoothT;
if (smoothT>1)
{
	for (i=1; i<=smoothT; i++)
	{
		for (j=1; j<=smoothT; j++)
		{
			if ((j==spinalPos)&&(i!=spinalPos)) { kernelBinT=kernelBinT+"1\t"; }
			if ((j==spinalPos)&&(i==spinalPos)) { kernelBinT=kernelBinT+centralPixel+"\t"; }
			if ((j!=spinalPos)&&(i!=spinalPos)) { kernelBinT=kernelBinT+"0\t"; }
			if ((j!=spinalPos)&&(i==spinalPos)) { kernelBinT=kernelBinT+"0\t"; }
		}
		kernelBinT=kernelBinT+"\n";
	}
	kernelBinT="text1=["+kernelBinT+"] stack";
}

//GET IMAGE STACK DETAILS
selectWindow("P");
getDimensions(width, height, channels, slices, frames);
timeGap=intTime/frames;

//DUPLICATE AND KEEP ORIGINALS FOR REFERENCE
selectWindow("P"); run("Duplicate...", "title=original_P duplicate");
selectWindow("S"); run("Duplicate...", "title=original_S duplicate");

//APPLYING MOVING AVERAGE OVER TIME AXIS ON ORIGINAL IMAGES
if (smoothT>1)
{
	for (iPS = 0; iPS < 2; iPS++)
	{
		//select a channel from "P" and "S"
		ChannelPS=substring("PS", iPS, iPS+1);
		selectWindow(ChannelPS); getPixelSize(DistanceUnit, pixelWidth, pixelHeight);
		selectWindow(ChannelPS); getMinAndMax(minPixelVal, maxPixelVal);
		selectWindow(ChannelPS); setVoxelSize(1, 1, 1, "pixel"); //removing scales for resclicing
		//CONVERT (X,Y,T) TO (X,T,Y)
		selectWindow(ChannelPS); run("Reslice [/]...", "output=1.000 start=Top avoid");
		selectWindow("Reslice of "+ChannelPS); rename(ChannelPS+"resliced");
		//APPLY MOVING AVERAGE OVER T USING PRE-DESIGNED KERNEL
		selectWindow(ChannelPS+"resliced"); run("Convolve...", kernelBinT);
		//CONVERT (X,T,Y) BACK TO (X,Y,T)
		selectWindow(ChannelPS+"resliced"); run("Reslice [/]...", "output=1.000 start=Top avoid");
		selectWindow("Reslice of "+ChannelPS+"resliced"); rename(ChannelPS+"reslicedBack");
		//GET RID OF OTHER IMAGES
		selectWindow(ChannelPS+"resliced"); close();
		selectWindow(ChannelPS); close();
		//RENAME AND RE-APPLY SCALES
		selectWindow(ChannelPS+"reslicedBack"); setVoxelSize(pixelWidth, pixelHeight, 1, DistanceUnit); setMinAndMax(minPixelVal, maxPixelVal);
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]"); run("Hyperstack to Stack");
		run("Duplicate...", "title=["+ChannelPS+"] duplicate range="+(1+(floor(smoothT/2)))+"-"+(frames-(floor(smoothT/2)))); //duplicating only the time frames with applied moving average
		selectWindow(ChannelPS+"reslicedBack"); close();
	}
}


//FIND THE FRAME WITH MAXIMUM INTENSITY IN EACH CHANNEL
selectWindow("P"); frames=nSlices;
maxP=0; maxS=0; maxPpos=0; maxSpos=0;
if (roiManager("count")>0) { roiManager("deselect"); roiManager("delete"); } //clear previously stored ROIs
selectWindow("P"); run("Select All"); roiManager("add"); //add full image count
run("Clear Results");
run("Set Measurements...", "mean integrated redirect=None decimal=9");
selectWindow("P"); roiManager("select", 0); roiManager("multi measure");
for (i=1; i<=frames; i++)
{
	count=parseFloat(getResultString("RawIntDen1", (i-1)));
	if (count>maxP) { maxP=count; maxPpos=i; }
}
run("Clear Results");
selectWindow("S"); roiManager("select", 0); roiManager("multi measure");
for (i=1; i<=frames; i++)
{
	count=parseFloat(getResultString("RawIntDen1", (i-1)));
	if (count>maxS) { maxS=count; maxSpos=i; }
}
roiManager("select", 0); roiManager("delete"); run("Clear Results");


//DUPLICATE THE IMAGE STACKS STARTING FROM MAX INTENSITY SLICE TO EQUAL FRAMES
startSlice=maxOf(maxPpos, maxSpos); numberOfslices=frames-startSlice;
selectWindow("P"); run("Duplicate...", "title=p duplicate range="+maxPpos+"-"+(maxPpos+numberOfslices));
selectWindow("S"); run("Duplicate...", "title=s duplicate range="+maxSpos+"-"+(maxSpos+numberOfslices));

//APPLY SMOOTHING ACCORDING TO USER INPUT
for (iPS = 0; iPS < 2; iPS++)
{
	//select a channel from "p" and "s"
	ChannelPS=substring("ps", iPS, iPS+1); selectWindow(ChannelPS);
	
	//Subtracting background
	run("Subtract...", "value="+background+" stack");
	
	//Applying pixel binning
	run("Convolve...", kernelBinXY);
}

selectWindow("s"); run("Multiply...", "value="+Gfactor+" stack"); //multiply Gfactor with the perpendicular channel
run("Duplicate...", "title=s2 duplicate");
selectWindow("s2");
run("Multiply...", "value=2");
//calculate numerator (I_p - Gfactor*I_s)
imageCalculator("Subtract create 32-bit stack", "p","s");
selectWindow("Result of p");
rename("num");
//calculate denominator (I_p + 2*Gfactor*I_s)
imageCalculator("Add create 32-bit stack", "p","s2");
selectWindow("Result of p");
rename("den");
//calculate r = (I_p - Gfactor*I_s)/(I_p + 2*Gfactor*I_s)
imageCalculator("Divide create 32-bit stack", "num","den");
selectWindow("Result of num");
rename("Anisotropy");

selectWindow("p"); close();
selectWindow("s"); close();
selectWindow("s2"); close();
selectWindow("num"); close();

selectWindow("P"); rename("P_smoothed");
selectWindow("S"); rename("S_smoothed");

if (keepSmooth==false)
{
	selectWindow("P_smoothed"); close();
	selectWindow("S_smoothed"); close();
}

close("ROI Manager"); close("Results");
selectWindow("original_P"); rename("P");
selectWindow("original_S"); rename("S");

selectWindow("den"); rename("Total Intensity"); run("Enhance Contrast", "saturated=0.35");
selectWindow("Anisotropy"); run("16 colors"); setMinAndMax(0.000000000, 0.400000000);