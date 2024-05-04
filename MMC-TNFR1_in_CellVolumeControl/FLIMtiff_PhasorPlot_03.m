% For Phasor diagram of FLIM Anisotropy data.
% Also Fluorescence Lifetime extraction from Phasor and fitting.
% uses .tif files as input.

clear all
%% Define experimental conditions and constants

nsTimeBin = 0.016; %Time Bin in NANOSECONDS

nHarmonic=1; %n-th order harmonics
decEndF=1e-3; %fraction of maximum total photon count that marks decay end
bgCountTF=0;%length of time to be considered as background at decay end
%set it 0 if not needed
%to manually define the end points of fluorescence emission decay
%decayStart=200; decayEnd=2500; %use manual input of start and end points

pxMASK=0; %intensity cutoff for object's pixels in total photon image
drawROI=1; %put =0 to not apply and use whole image, put =1 to apply

%Apply time-binning to reduce noise in the temporal decay of photon counts
tBIN=1; % add the photon count of tBIN bins
% tBin=1 means time binning is not going to be performed

%Apply 2D Kernel for pixel binning for spatial noise reduction
kSz=1; % implies (szk*szK) pixel binning (must be an odd number)
kWt=1; % weight multiplied to the central pixel (1 means pixels are summed)

%% Choose fitting model for Fluorescence Lifetime decay (optional)
%put =0 to not apply, put =1 to apply (fitting will slow the code)
fitE1=0;
fitE2=0;

%% Choose file to operate on
[fname,~] = uigetfile('*.tif','Select File');
if isequal(fname,0), disp('No file selected'); return, end

%% Read file
% Extract the main file name without the file type annotations
fnameMain=string(extractBetween(fname,1,(length(fname)-4)));

%read image stack
info = imfinfo(fname);%obtain all information about image stack
nTimeChannels = numel(info); %get number of images in stack
bD = info.BitDepth;
nMicronsPerPixel = info.XResolution; %resolution per image file
Simg1=imread(fname,1,'Info',info);%read 1 image slice
[nPixelsX,nPixelsY]=size(Simg1); %get image width/height
FLIMimg=zeros(nPixelsX,nPixelsY,nTimeChannels); %store the FLIM image
%loop to store whole image stack in a 3D matrix of DOUBLE type
for iimg=1:1:nTimeChannels
    Simg1=imread(fname,iimg,'Info',info);%read 1 image slice
    FLIMimg(:,:,iimg)=double(Simg1);%stored as DOUBLE matrix
end

res=nMicronsPerPixel/(1e-4); %image resolution in cm (Matlab default)

tic
%% Create intensity mask based on total photon count
MASKimg=sum(FLIMimg,3);
MASKimg(MASKimg<=pxMASK)=NaN;
MASKimg=(MASKimg.*0)+1;

%% Apply time binning
if (tBIN>1)
% Apply time binning
szTb=length(1:tBIN:nTimeChannels); %size of time binned stacks
ChPSb=zeros(nPixelsX,nPixelsY,szTb);
kT=1;
for iT=1:tBIN:(nTimeChannels-tBIN)
    ChPSb(:,:,kT)=sum(FLIMimg(:,:,iT:iT+tBIN),3);
    kT=kT+1;
end
%recalculate variables
nsTimeBin=nsTimeBin*tBIN; nTimeChannels=szTb;
FLIMimg=ChPSb;
clear ChPSb kT szTb
end

%% Apply pixel binning by convolution with 2D kernel
if (kSz>1)
% Define the 2D Kernel
cPxl=floor(kSz/2)+1; % central pixel position
Kernel=ones(kSz,kSz); %2D kernel (computes sum of szK*szK pixels)
Kernel(cPxl,cPxl)=Kernel(cPxl,cPxl)*kWt; %adding weighted central pixel
% Apply pixel binning to parallel and perpendicular channels
for iK=1:nTimeChannels
    FLIMimg(:,:,iK)=conv2(FLIMimg(:,:,iK),Kernel,'same');
end
end

%% Determine the time point signifying the end of the decay curve (FLIM)
if ((exist('decayStart','var')==1)&&(exist('decayEnd','var')==1))
    %using manual input
    maxIpos=decayStart; minIpos=decayEnd;
else
    % Find the maximum and minimum Intensity from total photon emission
    SumI=sum(FLIMimg,[1,2]); SumI=reshape(SumI,[1,nTimeChannels]);
    SumI=smoothdata(SumI,'sgolay');
    [maxSumI,maxSumIpos]=max(SumI);
    SumIdecay=SumI(maxSumIpos:nTimeChannels);
    [~,minDecayPoint]=min(abs(SumIdecay-(decEndF*maxSumI)));
    minDecayPoint=minDecayPoint+maxSumIpos-1;
    maxIpos=maxSumIpos;
    minIpos=minDecayPoint;
    clear SumI SumIdecay maxSumI maxSumIpos minDecayPoint
end

szT=(minIpos-maxIpos)+1; %size of time frames to be considered for decay

%% Separate the Decay part from the total FLIM image stack
DecayPart=FLIMimg(:,:,maxIpos:minIpos); %Total intensity decay (for FLIM)
TotalPhotons=sum(DecayPart,3);

%% Substract background counts at the end of the decay (if needed)
if (bgCountTF>0)
bgHalf=floor(bgCountTF/2);
%define and subtract the mean photon counts of background time frames
bgPS=mean(FLIMimg(:,:,(minIpos-bgHalf):(minIpos+bgHalf)),3);
DecayPart=DecayPart-bgPS;
clear bgPS
end

%% Find region of interest interactively using the total photon image
if (drawROI==1)
TotPhotMax=max(max(TotalPhotons)); TotPhotMin=min(min(TotalPhotons));
PhotonSumImg=(TotalPhotons-TotPhotMin)*65535/(TotPhotMax-TotPhotMin);
ROIfig=figure; imshow(uint16(PhotonSumImg)); colormap('turbo');
%
xlabel('Draw freehand ROI denoting Object');
ROIobj=drawfreehand('LineWidth',1,'Color','white');
%
%{
xlabel('Draw Polygonal ROI denoting Object');
ROIobj=drawpolygon('LineWidth',1,'Color','white');
%}
ROImask=createMask(ROIobj);
ROImask=double(ROImask);
ROImask(ROImask==0)=NaN; %Convert unselected regions to NaN
%Apply user drawn mask on relevant image stacks
DecayPart=bsxfun(@times,DecayPart,cast(ROImask,class(DecayPart)));
close(ROIfig);
clear ROIfig ROIobj
end

%% Apply Intensity mask on whole image stacks
DecayPart=bsxfun(@times,DecayPart,cast(MASKimg,class(DecayPart)));

%% Phasor calculation
% Calculate the angular frequency (omega) for Phasor calculation
timePeriod=(minIpos-maxIpos)*nsTimeBin*(10^(-9)); %in s
omega=2*pi*nHarmonic/timePeriod; %in Hz
% Create Time, Cosine, and Sine vectors
Tns=0:nsTimeBin:(nsTimeBin*(szT-1)); %Time vector in ns
Ts=Tns*(10^(-9)); %Time vector in s
Vcos=double(cos(omega*Ts));
Vsin=double(sin(omega*Ts));
% Define Cosine and Sine matrices
Mcos=ones(nPixelsX,nPixelsY,szT);
Msin=ones(nPixelsX,nPixelsY,szT);
% Multiply Vcos, Vsin to the time axes of Mcos, Msin at every pixel
Mcos=bsxfun(@times,Mcos,reshape(Vcos,1,1,szT));
Msin=bsxfun(@times,Msin,reshape(Vsin,1,1,szT));
% Calculate Integrations (Summations)
MG=DecayPart.*Mcos;
MS=DecayPart.*Msin;
% Calculate Phasors and Fluorescence Lifetime
G=sum(MG,3)./TotalPhotons;
S=sum(MS,3)./TotalPhotons;
%Linearize G and S for plotting later
LinearG=reshape(G,[],1); LinearG=LinearG(~isnan(LinearG));
LinearS=reshape(S,[],1); LinearS=LinearS(~isnan(LinearS));
TauP=((10^9)/omega)*(S./G); %in ns
% Find magnitude and phase angle of the phasor points (vectors)
PhMag=sqrt((G.^2)+(S.^2));
PhAng=atan2(S,G);

%% Plotting histogram of the Phasor points on the universal semi-circle
%figure,plot(G,S,'k.','MarkerSize',1); xlim([-0.1 1.1]); ylim([0 0.7]);
%set(gca,'XColor','none','YColor','none'); set(gca,'color','none');
figure('units','normalized','position',[0.0 0.2 1.0 0.45])
%plot the universal semi-circle for reference
Gusc=0:0.0001:1; Susc=sqrt(Gusc-Gusc.^2); %define the universal semi-circle
grayline=[0.7 0.7 0.7]; %line color of the plotted universal semi-circle
subplot(1,2,1)
plot(Gusc,Susc,'-','MarkerSize',2);
xlabel('G \rightarrow'); ylabel('S \rightarrow');
colororder(grayline)
hold on
%plot a lifetime scale bar (in ns) on the universal semi-circle
TauScaleUSC=Susc./(Gusc*omega*(10^-9)); %Lifetime scale (in ns)
nPhasorMarkers=20;
%cTS=jet(nPoints);
for iTSB=1:1:nPhasorMarkers
    TickPos=TauScaleUSC-iTSB;
    [~,mTPp]=min(abs(TickPos));
    plPh=plot(Gusc(mTPp),Susc(mTPp),'o');
    plPh.MarkerSize=8;
    plPh.MarkerFaceColor='w'; %cTS(iTSB,:);
    plPh.MarkerEdgeColor='k';
    labelText=sprintf('%d',iTSB);
    markerLabel=text(Gusc(mTPp),Susc(mTPp),labelText);
    markerLabel.Color='k'; %1-cTS(iTSB,:);
    markerLabel.FontSize=(plPh.MarkerSize-2);
    markerLabel.HorizontalAlignment='center';
    markerLabel.VerticalAlignment='middle';
    hold on
end
hold on
% Plot histogram of Phasors on the universal semi-circle
nPxB=1; %for  n_P*n_P  pixel binning of phasor points
nbins=round(nPixelsX/nPxB);
Hist2Dplot=histogram2(G,S,nbins); %create bivariate histogram
Hist2Dplot.DisplayStyle='tile'; %show as tiles, not bars in 3D plot
Hist2Dplot.ShowEmptyBins='off'; %don't show empty bins
Hist2Dplot.EdgeColor='none'; %don't show edge colors for a smooth view
view(2) %view in 2D
shading interp %interpolated shading
cbarPhasor=colorbar; %show colorbar
cbarPhasor.Label.String='Number of Pixels';
title('Phasor Plot Histogram');
xlim([-0.1 1.1]); ylim([0 0.7]);
hold off
% Plot contour map of Phasors on the universal semi-circle
%
% Define bin edges and bins for (G,S) values (the Phasors)
%Gedges=linspace(-0.1,1.1,nBinsC); Sedges=linspace(-0.1,1.1,nBinsC);
limG1=min(min(G)); limG2=max(max(G));
limS1=min(min(S)); limS2=max(max(S));
Gedges=linspace(limG1,limG2,nbins);
Sedges=linspace(limS1,limS2,nbins);
%obtain 2D histogram
HistGS=histcounts2(G,S,Gedges,Sedges);
%use the bin edges to find the bin centers in the G and S plane
[GG,SS]=meshgrid(mean([Gedges(1:end-1);Gedges(2:end)],1),...
    mean([Sedges(1:end-1);Sedges(2:end)],1));
%plot 2D histogram with contourlines
subplot(1,2,2)
%ContourGS=contourf(GG,SS,HistGS','LineColor','none');
ContourGS=contour(GG,SS,HistGS');
cmap=turbo;
cmap(1,:)=[1,1,1]; %white background
%cmap(1,:)=[0,0,0]; %black background
colormap(cmap);
shading interp %interpolated shading
cbarPhasor=colorbar; %show colorbar
cbarPhasor.Label.String='Number of Pixels';
title('Zoomed in Contour Map of Phasor Points');
xlabel('G \rightarrow'); ylabel('S \rightarrow');
%

%% Fit photon emission decay for each individual pixel
%to store tau values from fitting
if (fitE1>0), TauE1=zeros(nPixelsX,nPixelsY); end
if (fitE2>0), TauE2=zeros(nPixelsX,nPixelsY); end

if ((fitE1==1)||(fitE2==1))
for x=1:nPixelsX
    parfor y=1:nPixelsY
        %extract TCSPC I(t) of the pixel at (x,y)
        I=reshape((double(DecayPart(x,y,:))),[1,szT]);
        nanTEST=isnan(sum(I));
        % FIT Intensity Decay BY Mono-Exponential
        if (fitE1~=0)
        if (nanTEST~=1)
            %fit to mono-exponential decay
            VarI1=[1,-1/TauP(x,y)];
            [e1fit,gof]=fit(Tns',I','exp1','StartPoint',VarI1);
            %get Coefficients
            Cf1=coeffvalues(e1fit);
            TauE1(x,y)=abs(1/Cf1(2));
            %extract tau only if the fitting is good
            %{
            if (gof.adjrsquare>=0.95)
                TauE1(x,y)=abs(1/Cf1(2));
            else
                TauE1(x,y)=NaN;
            end
            %}
        else
            TauE1(x,y)=NaN;
        end
        end
        % FIT Intensity Decay BY Bi-Exponential
        if (fitE2~=0)
        if (nanTEST~=1)
            %fit to bi-exponential decay
            VarI2=[1,-1/TauP(x,y),1,-1/TauP(x,y)];
            [e2fit,gof]=fit(Tns',I','exp2','StartPoint',VarI2);
            %get Coefficients
            Cf2=coeffvalues(e2fit);
            %{
            %INTENSITY-AVERAGED LIFETIME
            num=(Cf2(1)*(1/Cf2(2))^2)+(Cf2(3)*(1/Cf2(4))^2);
            den=(Cf2(1)*abs(1/Cf2(2)))+(Cf2(3)*abs(1/Cf2(4)));
            %}
            %
            %AMPLITUDE-AVERAGED LIFETIME
            num=(Cf2(1)*abs(1/Cf2(2)))+(Cf2(3)*abs(1/Cf2(4)));
            den=Cf2(1)+Cf2(3);
            %
            TauE2(x,y)=num/den;
            %extract tau only if the fitting is good
            %{
            if (gof.adjrsquare>=0.95)
                %{
                %INTENSITY-AVERAGED LIFETIME
                num=(Cf2(1)*(1/Cf2(2))^2)+(Cf2(3)*(1/Cf2(4))^2);
                den=(Cf2(1)*abs(1/Cf2(2)))+(Cf2(3)*abs(1/Cf2(4)));
                %}
                %
                %AMPLITUDE-AVERAGED LIFETIME
                num=(Cf2(1)*abs(1/Cf2(2)))+(Cf2(3)*abs(1/Cf2(4)));
                den=Cf2(1)+Cf2(3);
                %
            else
                TauE2(x,y)=NaN;
            end
            %}
        else
            TauE2(x,y)=NaN;
        end
        end
    end
end
end

%% Show the fluorescence lifetime map
%{
figure('units','normalized','position',[0.25 0.1 0.2 0.6])
imshow(TauP)
title('Lifetime map (Phasor Plot)');
cmap=turbo;
cmap(1,:)=[1,1,1]; %white background
%cmap(1,:)=[0,0,0]; %black background
colormap(cmap);
caxis([min(min(TauP)) max(max(TauP))])
cbarTau=colorbar;
cbarTau.Label.String='Lifetime (ns)';
%}

%% SAVE RESULTS AS TIFF IMAGES
%
%save Intensity Map as a 32bit .tif using TIFF class
sfname=strcat('TotalIntensityMap-',fnameMain,'.tif');
Tiff_IntMap = Tiff(char(sfname),'w');
setTag(Tiff_IntMap,'Photometric',Tiff.Photometric.MinIsBlack);
setTag(Tiff_IntMap,'Compression',Tiff.Compression.None);
setTag(Tiff_IntMap,'BitsPerSample',32);
setTag(Tiff_IntMap,'SamplesPerPixel',1);
setTag(Tiff_IntMap,'SampleFormat',Tiff.SampleFormat.IEEEFP);
setTag(Tiff_IntMap,'ImageLength',nPixelsX);
setTag(Tiff_IntMap,'ImageWidth',nPixelsY);
setTag(Tiff_IntMap,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
setTag(Tiff_IntMap,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
setTag(Tiff_IntMap,'XResolution',res);
setTag(Tiff_IntMap,'YResolution',res);
write(Tiff_IntMap,single(TotalPhotons));
Tiff_IntMap.close();

%save FLIM Map (Phasor) as a 32bit .tif using TIFF class
sfname=strcat('TauMap-byPhasor-',fnameMain,'.tif');
Tiff_TauPh = Tiff(char(sfname),'w');
setTag(Tiff_TauPh,'Photometric',Tiff.Photometric.MinIsBlack);
setTag(Tiff_TauPh,'Compression',Tiff.Compression.None);
setTag(Tiff_TauPh,'BitsPerSample',32);
setTag(Tiff_TauPh,'SamplesPerPixel',1);
setTag(Tiff_TauPh,'SampleFormat',Tiff.SampleFormat.IEEEFP);
setTag(Tiff_TauPh,'ImageLength',nPixelsX);
setTag(Tiff_TauPh,'ImageWidth',nPixelsY);
setTag(Tiff_TauPh,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
setTag(Tiff_TauPh,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
setTag(Tiff_TauPh,'XResolution',res);
setTag(Tiff_TauPh,'YResolution',res);
write(Tiff_TauPh,single(TauP));
Tiff_TauPh.close();

%save FLIM Map (from fitting) as a 32bit .tif using TIFF class
if (fitE1>0)
sfname=strcat('TauMap-byExp1fit-',fnameMain,'.tif');
Tiff_TauExp1 = Tiff(char(sfname),'w');
setTag(Tiff_TauExp1,'Photometric',Tiff.Photometric.MinIsBlack);
setTag(Tiff_TauExp1,'Compression',Tiff.Compression.None);
setTag(Tiff_TauExp1,'BitsPerSample',32);
setTag(Tiff_TauExp1,'SamplesPerPixel',1);
setTag(Tiff_TauExp1,'SampleFormat',Tiff.SampleFormat.IEEEFP);
setTag(Tiff_TauExp1,'ImageLength',nPixelsX);
setTag(Tiff_TauExp1,'ImageWidth',nPixelsY);
setTag(Tiff_TauExp1,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
setTag(Tiff_TauExp1,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
setTag(Tiff_TauExp1,'XResolution',res);
setTag(Tiff_TauExp1,'YResolution',res);
write(Tiff_TauExp1,single(TauE1));
Tiff_TauExp1.close();
end

%save FLIM Map (from fitting) as a 32bit .tif using TIFF class
if (fitE2>0)
sfname=strcat('TauMap-byExp2fit-',fnameMain,'.tif');
Tiff_TauExp2 = Tiff(char(sfname),'w');
setTag(Tiff_TauExp2,'Photometric',Tiff.Photometric.MinIsBlack);
setTag(Tiff_TauExp2,'Compression',Tiff.Compression.None);
setTag(Tiff_TauExp2,'BitsPerSample',32);
setTag(Tiff_TauExp2,'SamplesPerPixel',1);
setTag(Tiff_TauExp2,'SampleFormat',Tiff.SampleFormat.IEEEFP);
setTag(Tiff_TauExp2,'ImageLength',nPixelsX);
setTag(Tiff_TauExp2,'ImageWidth',nPixelsY);
setTag(Tiff_TauExp2,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
setTag(Tiff_TauExp2,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
setTag(Tiff_TauExp2,'XResolution',res);
setTag(Tiff_TauExp2,'YResolution',res);
write(Tiff_TauExp2,single(TauE2));
Tiff_TauExp2.close();
end

%save the whole FLIM image stack as a 32bit .tif using TIFF class
%{
sfname=strcat('FLIMstack-',fnameMain,'.tif');
Tiff_FLIM = Tiff(char(sfname),'w');
setTag(Tiff_FLIM,'Photometric',Tiff.Photometric.MinIsBlack);
setTag(Tiff_FLIM,'Compression',Tiff.Compression.None);
setTag(Tiff_FLIM,'BitsPerSample',32);
setTag(Tiff_FLIM,'SamplesPerPixel',nTimeChannels);
setTag(Tiff_FLIM,'SampleFormat',Tiff.SampleFormat.IEEEFP);
setTag(Tiff_FLIM,'ImageLength',nPixelsX);
setTag(Tiff_FLIM,'ImageWidth',nPixelsY);
setTag(Tiff_FLIM,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
setTag(Tiff_FLIM,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
setTag(Tiff_FLIM,'XResolution',res);
setTag(Tiff_FLIM,'YResolution',res);
write(Tiff_FLIM,single(ChPS));
Tiff_FLIM.close();
%}

%}
toc