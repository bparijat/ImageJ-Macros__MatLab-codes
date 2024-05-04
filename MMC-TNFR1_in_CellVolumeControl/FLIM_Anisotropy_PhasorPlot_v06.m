% For Phasor diagram of FLIM Anisotropy data.
% Also Fluorescence Lifetime extraction from Phasor and fitting.
% uses .bin files as input.
% The PARALLEL and PERPENDICULAR channels should have similar file names..
% ..but "P" (for Parallel) or "S" (for Perpendicular) at the end.
% eg: Sample1_c1P.bin, Sample1_c2S.bin or Sample1_c2P.bin, Sample1_c1S.bin
% Both files to be operated on must be in the same folder.

clear all
%% Define experimental conditions and constants
nHarmonic=1; %n-th order harmonics
decEndF=1e-3; %fraction of maximum total photon count that marks decay end
bgCountTF=50;%length of time to be considered as background at decay end
%manually define the end points of fluorescence emission decay
%decayStart=200; decayEnd=2500; %use manual input of start and end points

gFactor=3; %G-factor, multiply with the Perpendicular channel (ChS)
%gFactor=0.43; %G-factor, multiply with the Perpendicular channel (ChS)
fitAni=0; %put =1 to fit time resolved anisotropy decay 
nFani=300; %number of time frames for fitting anisotropy decay

pxMASK=0; %intensity cutoff for object's pixels in total photon image
drawROI=0; %put =0 to not apply, put =1 to apply

%Apply time-binning
tBIN=1; % add the photon count of tBIN time frames to make one frame

%Apply 2D Kernel for pixel binning
kSz=1; % implies (szk*szK) pixel binning (must be an odd number)
kWt=1; % weight multiplied to the central pixel (1 means pixels are summed)

%% Choose fitting model for Fluorescence Lifetime decay
%put =0 to not apply, put =1 to apply
fitE1=0;
fitE2=0;

%
%% Choose filenames of the parallel and perpendicular channels
[fnameP,~] = uigetfile('*P.bin','Select Parallel Channel');
if isequal(fnameP,0), disp('No file selected'); return, end
[fnameS,~] = uigetfile('*S.bin','Select Perpendicular Channel');
if isequal(fnameS,0), disp('No file selected'); return, end

%% Read files from folder for operation
% Extract the main file name without the file type and channel annotations
fnameMain=string(extractBetween(fnameP,1,(length(fnameP)-8)));

% Open and read headers (image size, time channels etc)
%has to be in this exact sequence
fP = fopen(fnameP); %reading Parallel channel
nPixelsXP = fread(fP, 1, 'int32');
nPixelsYP = fread(fP, 1, 'int32');
nMicronsPerPixelP = fread(fP, 1, 'float32');
nTimeChannelsP = fread(fP, 1, 'int32');
nsTimeBinP = fread(fP, 1, 'float32');
%create the Parallel channel image stack
ChP=zeros(nPixelsXP,nPixelsYP,nTimeChannelsP);
for x=1:1:nPixelsXP
    for y=1:1:nPixelsYP
        ChP(x,y,:)=fread(fP, nTimeChannelsP, 'int32');
    end
end
fclose(fP);
fS = fopen(fnameS); %reading Perpendicular channel
nPixelsXS = fread(fS, 1, 'int32');
nPixelsYS = fread(fS, 1, 'int32');
nMicronsPerPixelS = fread(fS, 1, 'float32');
nTimeChannelsS = fread(fS, 1, 'int32');
nsTimeBinS = fread(fS, 1, 'float32');
%create the Perpendicular channel image stack
ChS=zeros(nPixelsXS,nPixelsYS,nTimeChannelsS);
for x=1:1:nPixelsXS
    for y=1:1:nPixelsYS
        ChS(x,y,:)=fread(fS, nTimeChannelsS, 'int32');
    end
end
fclose(fS);

% Check if imaging parameters are equal for both channels
if (nPixelsXP==nPixelsXS), nPixelsX=nPixelsXP; else, return; end
if (nPixelsYP==nPixelsYS), nPixelsY=nPixelsYP; else, return; end
if (nMicronsPerPixelP==nMicronsPerPixelS)
    nMicronsPerPixel=nMicronsPerPixelP;
else
    return;
end
if (nTimeChannelsP==nTimeChannelsS)
    nTimeChannels=nTimeChannelsP;
else
    return;
end
if (nsTimeBinP==nsTimeBinS)
    nsTimeBin=nsTimeBinP;
else
    return;
end
clear nPixelsXP nPixelsYP nMicronsPerPixelP nTimeChannelsP nsTimeBinP
clear nPixelsXS nPixelsYS nMicronsPerPixelS nTimeChannelsS nsTimeBinS

res=(1/(nMicronsPerPixel*1e-4)); %image resolution in cm

%% Create the total intensity FLIM decay image
ChPS=double(ChP)+double(ChS);
%

tic
%% Create intensity mask based on total photon count
MASKimg=sum(ChPS,3);
MASKimg(MASKimg<=pxMASK)=NaN;
MASKimg=(MASKimg.*0)+1;

%% Apply time binning
if (tBIN>1)
% Apply time binning
szTb=length(1:tBIN:nTimeChannels); %size of time binned stacks
ChPSb=zeros(nPixelsX,nPixelsY,szTb);
ChPb=zeros(nPixelsX,nPixelsY,szTb);
ChSb=zeros(nPixelsX,nPixelsY,szTb);
kT=1;
for iT=1:tBIN:(nTimeChannels-tBIN)
    ChPSb(:,:,kT)=sum(ChPS(:,:,iT:iT+tBIN),3);
    ChPb(:,:,kT)=sum(ChP(:,:,iT:iT+tBIN),3);
    ChSb(:,:,kT)=sum(ChS(:,:,iT:iT+tBIN),3);
    kT=kT+1;
end
%recalculate variables
nsTimeBin=nsTimeBin*tBIN; nTimeChannels=szTb;
ChPS=ChPSb; ChP=ChPb; ChS=ChSb;
clear ChPSb ChPb ChSb kT szTb
end

%% Apply pixel binning by convolution with 2D kernel
if (kSz>1)
% Define the 2D Kernel
cPxl=floor(kSz/2)+1; % central pixel position
Kernel=ones(kSz,kSz); %2D kernel (computes sum of szK*szK pixels)
Kernel(cPxl,cPxl)=Kernel(cPxl,cPxl)*kWt; %adding weighted central pixel
% Apply pixel binning to parallel and perpendicular channels
for iK=1:nTimeChannels
    ChPS(:,:,iK)=conv2(ChPS(:,:,iK),Kernel,'same');
    ChP(:,:,iK)=conv2(ChP(:,:,iK),Kernel,'same');
    ChS(:,:,iK)=conv2(ChS(:,:,iK),Kernel,'same');
end
end

%% Determine the time point signifying the end of the decay curve (FLIM)
if ((exist('decayStart','var')==1)&&(exist('decayEnd','var')==1))
    %using manual input
    maxIpos=decayStart; minIpos=decayEnd;
else
    % Find the maximum and minimum Intensity from total photon emission
    SumI=sum(ChPS,[1,2]); SumI=reshape(SumI,[1,nTimeChannels]);
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

%% Find individual maxima position for parallel and perpendicular channels
SumP=sum(ChP,[1,2]); SumP=reshape(SumP,[1,nTimeChannels]);
SumP=smoothdata(SumP,'sgolay');
[~,maxSumPpos]=max(SumP); %maxima position of parallel channel
SumS=sum(ChS,[1,2]); SumS=reshape(SumS,[1,nTimeChannels]);
SumS=smoothdata(SumS,'sgolay');
[~,maxSumSpos]=max(SumS); %maxima position of perpendicular channel
clear SumP SumS

%% Separate the Decay part from the total FLIM image stack
DecayPS=ChPS(:,:,maxIpos:minIpos); %Total intensity decay (for FLIM)
DecayP=ChP(:,:,maxSumPpos:(maxSumPpos+szT-1)); %Parallel channel decay
DecayS=ChS(:,:,maxSumSpos:(maxSumSpos+szT-1)); %Perpendicular channel decay
TotalPhotons=sum(DecayPS,3);

%% Substract background counts (found at the end of the decay)
if (bgCountTF>0)
bgHalf=floor(bgCountTF/2);
%define and subtract the mean photon counts of background time frames
%for parallel+perpendicular channel
bgPS=mean(ChPS(:,:,(minIpos-bgHalf):(minIpos+bgHalf)),3);
DecayPS=DecayPS-bgPS;
% for parallel channel
bgP=mean(ChP(:,:,(maxSumPpos+szT-1-bgHalf):(maxSumPpos+szT-1+bgHalf)),3);
DecayP=DecayP-bgP;
% for perpendicular channel
bgS=mean(ChS(:,:,(maxSumSpos+szT-1-bgHalf):(maxSumSpos+szT-1+bgHalf)),3);
DecayS=DecayS-bgS;
clear bgPS bgP bgS
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
DecayPS=bsxfun(@times,DecayPS,cast(ROImask,class(DecayPS)));
DecayP=bsxfun(@times,DecayP,cast(ROImask,class(DecayP)));
DecayS=bsxfun(@times,DecayS,cast(ROImask,class(DecayS)));
close(ROIfig);
clear ROIfig ROIobj
end

%% Apply Intensity mask on whole image stacks
DecayPS=bsxfun(@times,DecayPS,cast(MASKimg,class(DecayPS)));
DecayP=bsxfun(@times,DecayP,cast(MASKimg,class(DecayP)));
DecayS=bsxfun(@times,DecayS,cast(MASKimg,class(DecayS)));

%% Generate steady state fluorescence anisotropy image
IP=sum(DecayP,3); IS=sum(DecayS,3);
AniSS=(IP-gFactor*IS)./(IP+2*gFactor*IS);

%% Generate time resolved fluorescence anisotropy image stack from raw data
AniTR=(DecayP-gFactor*DecayS)./(DecayP+2*gFactor*DecayS);

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
MG=DecayPS.*Mcos;
MS=DecayPS.*Msin;
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
figure('units','normalized','position',[0.0 0.1 0.25 0.6])
%plot the universal semi-circle for reference
Gusc=0:0.0001:1; Susc=sqrt(Gusc-Gusc.^2); %define the universal semi-circle
grayline=[0.7 0.7 0.7]; %line color of the plotted universal semi-circle
subplot(2,1,1)
plot(Gusc,Susc,'-','MarkerSize',2);
xlabel('G \rightarrow'); ylabel('S \rightarrow');
colororder(grayline)
hold on
%plot a lifetime scale bar (in ns) on the universal semi-circle
TauScaleUSC=Susc./(Gusc*omega*(10^-9)); %Lifetime scale (in ns)
nPhasorMarkers=10;
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
subplot(2,1,2)
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
Theta=zeros(nPixelsX,nPixelsY); %rotational correlation time
AniLim=zeros(nPixelsX,nPixelsY); %limiting anisotropy (r0)
Xr=Tns(1:nFani);
%to store tau values from fitting
if (fitE1>0), TauE1=zeros(nPixelsX,nPixelsY); end
if (fitE2>0), TauE2=zeros(nPixelsX,nPixelsY); end

if ((fitAni==1)||(fitE1==1)||(fitE2==1))
for x=1:nPixelsX
    parfor y=1:nPixelsY
        %% Fit fluorescence anisotropy decay for each individual pixel
        %extract anisotropy decay of the pixel at (x,y)
        A=reshape((double(AniTR(x,y,:))),[1,szT]);
        nanTEST=isnan(sum(A));
        if (nanTEST~=1)
            % Apply smoothing to decay curves
            %A=smoothdata(A,'sgolay');
            A=medfilt1(A);
            % Select decay section to fit the curve
            Yr=A(1:nFani);
            % Initialize variable start points for fitting
            FitStart=[0.4,(1/20)];
            % Fit to mono-exponential decay
            [rfit,gof]=fit(Xr',Yr','exp1','StartPoint',FitStart);
            %get Coefficients
            Cfr=coeffvalues(rfit);
            Theta(x,y)=abs(1/Cfr(2));
            AniLim(x,y)=Cfr(1);
        else
            Theta(x,y)=NaN;
            AniLim(x,y)=NaN;
        end
        
        %% Fit Intensity Decay for each individual pixel
        %extract TCSPC I(t) of the pixel at (x,y)
        I=reshape((double(DecayPS(x,y,:))),[1,szT]);
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

%% Show the fluorescence lifetime and anisotropy maps
figure('units','normalized','position',[0.25 0.1 0.2 0.6])
subplot(2,1,1)
imshow(TauP)
title('Lifetime map (Phasor Plot)');
cmap=turbo;
cmap(1,:)=[1,1,1]; %white background
%cmap(1,:)=[0,0,0]; %black background
colormap(cmap);
caxis([min(min(TauP)) max(max(TauP))])
cbarTau=colorbar;
cbarTau.Label.String='Lifetime (ns)';
subplot(2,1,2)
imshow(AniSS)
title('Steady State Anisotropy map');
cmap=turbo;
cmap(1,:)=[1,1,1]; %white background
%cmap(1,:)=[0,0,0]; %black background
colormap(cmap);
caxis([0.2,0.4])
cbarTau=colorbar;
cbarTau.Label.String='Anisotropy';

if (fitAni==1)
figure('units','normalized','position',[0.45 0.1 0.2 0.6])
subplot(2,1,1)
imshow(AniLim)
title('r_0 map');
cmap=turbo;
cmap(1,:)=[1,1,1]; %white background
%cmap(1,:)=[0,0,0]; %black background
colormap(cmap);
%caxis([min(min(AniLim)),0.4])
caxis([0.2,0.4])
cbarTau=colorbar;
cbarTau.Label.String='r_0';
subplot(2,1,2)
imshow(Theta)
title('\theta map');
cmap=turbo;
cmap(1,:)=[1,1,1]; %white background
%cmap(1,:)=[0,0,0]; %black background
colormap(cmap);
%caxis([min(min(Theta)),50])
caxis([0,50])
cbarTau=colorbar;
cbarTau.Label.String='Rotational Correlation time (ns)';
end

%% Test the nature of the decay curves
if (fitAni==1)
Zr=nanmean(AniTR,[1,2]);
Zr=reshape(Zr,[1,szT]); Zr=Zr';%define the anisotropy decay as y-column
Zs=medfilt1(Zr); %smooth data by median filter
Zx=1:1:szT; Zx=Zx'; %define time points as x-column
figure('units','normalized','position',[0.65 0.1 0.2 0.6])
subplot(2,1,1)
plot(Zx,Zr,Zx,Zs), xlim([0,Zx(end)]), ylim([-2,2]);
title('Mean Anisotropy Decay of all pixels');
legend('Raw','Smoothed','Location','northwest');
subplot(2,1,2)
plot(Zx,Zr,Zx,Zs), xlim([0,nFani]), ylim([-0.05,0.6]);
title('Anisotropy Decay section used for fitting');
legend('Raw','Smoothed','Location','northwest');
end

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

%save the Steady State Anisotropy map as a 32bit .tif
sfname=strcat('AnisotropySS-',fnameMain,'.tif');
Tiff_AniSS = Tiff(char(sfname),'w');
setTag(Tiff_AniSS,'Photometric',Tiff.Photometric.MinIsBlack);
setTag(Tiff_AniSS,'Compression',Tiff.Compression.None);
setTag(Tiff_AniSS,'BitsPerSample',32);
setTag(Tiff_AniSS,'SamplesPerPixel',1);
setTag(Tiff_AniSS,'SampleFormat',Tiff.SampleFormat.IEEEFP);
setTag(Tiff_AniSS,'ImageLength',nPixelsX);
setTag(Tiff_AniSS,'ImageWidth',nPixelsY);
setTag(Tiff_AniSS,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
setTag(Tiff_AniSS,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
setTag(Tiff_AniSS,'XResolution',res);
setTag(Tiff_AniSS,'YResolution',res);
write(Tiff_AniSS,single(AniSS));
Tiff_AniSS.close();

if (fitAni==1)
%save the Limiting Anisotropy (r0) map as a 32bit .tif
sfname=strcat('r0map-',fnameMain,'.tif');
Tiff_AniLim = Tiff(char(sfname),'w');
setTag(Tiff_AniLim,'Photometric',Tiff.Photometric.MinIsBlack);
setTag(Tiff_AniLim,'Compression',Tiff.Compression.None);
setTag(Tiff_AniLim,'BitsPerSample',32);
setTag(Tiff_AniLim,'SamplesPerPixel',1);
setTag(Tiff_AniLim,'SampleFormat',Tiff.SampleFormat.IEEEFP);
setTag(Tiff_AniLim,'ImageLength',nPixelsX);
setTag(Tiff_AniLim,'ImageWidth',nPixelsY);
setTag(Tiff_AniLim,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
setTag(Tiff_AniLim,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
setTag(Tiff_AniLim,'XResolution',res);
setTag(Tiff_AniLim,'YResolution',res);
write(Tiff_AniLim,single(AniLim));
Tiff_AniLim.close();
%save the rotational correlation time map as a 32bit .tif
sfname=strcat('ThetaMap-ns-',fnameMain,'.tif');
Tiff_Theta = Tiff(char(sfname),'w');
setTag(Tiff_Theta,'Photometric',Tiff.Photometric.MinIsBlack);
setTag(Tiff_Theta,'Compression',Tiff.Compression.None);
setTag(Tiff_Theta,'BitsPerSample',32);
setTag(Tiff_Theta,'SamplesPerPixel',1);
setTag(Tiff_Theta,'SampleFormat',Tiff.SampleFormat.IEEEFP);
setTag(Tiff_Theta,'ImageLength',nPixelsX);
setTag(Tiff_Theta,'ImageWidth',nPixelsY);
setTag(Tiff_Theta,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
setTag(Tiff_Theta,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
setTag(Tiff_Theta,'XResolution',res);
setTag(Tiff_Theta,'YResolution',res);
write(Tiff_Theta,single(Theta));
Tiff_Theta.close();
end

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

%save the Anisotropy Decay part as a 32bit .tif image stack
%{
sfname=strcat('AnisotropyTR-',fnameMain,'.tif');
Tiff_AniTR = Tiff(char(sfname),'w');
setTag(Tiff_AniTR,'Photometric',Tiff.Photometric.MinIsBlack);
setTag(Tiff_AniTR,'Compression',Tiff.Compression.None);
setTag(Tiff_AniTR,'BitsPerSample',32);
setTag(Tiff_AniTR,'SamplesPerPixel',szT);
setTag(Tiff_AniTR,'SampleFormat',Tiff.SampleFormat.IEEEFP);
setTag(Tiff_AniTR,'ImageLength',nPixelsX);
setTag(Tiff_AniTR,'ImageWidth',nPixelsY);
setTag(Tiff_AniTR,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
setTag(Tiff_AniTR,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
setTag(Tiff_AniTR,'XResolution',res);
setTag(Tiff_AniTR,'YResolution',res);
write(Tiff_AniTR,single(AniTR));
Tiff_AniTR.close();
%}

%}
toc