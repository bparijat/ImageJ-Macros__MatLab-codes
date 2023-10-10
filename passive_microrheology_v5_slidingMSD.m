clear all
bR=1.0; bG=1.0; bB=1.0; %custom background colour for plotting

%% read all files from folder
FileNames=dir('Re*.xls');  %array stores file names in current folder
nf=length(FileNames);   %stores number of files read
TjN=zeros(nf,1);

%% do or don't (1=perform, 0=don't perform) [don't perform to save time]
displayPlots=1; %dispay all graphs for mean MSD, D, G',G" etc per file
plotAllTraj=0; %plot all trajectories
subtractRefTraj=0;  %subtract first trajectory to eliminate drift
ref1=1; %reference trajectory for subtraction (trajectory array index)
plotSubtractedTraj=1; %plot subtracted trajectories
if (subtractRefTraj==0), plotSubtractedTraj=0; ref1=0; end
saveMeanRheologyCurves=1; %save mean MSD, D, G',G" etc for each file
saveTwoFrequencyValues=1; %save quantities at 10Hz, 1Hz etc
logPlot=1; slogPlot=0; %choose between EITHER log-log plot or semi-log plot
if (logPlot==1), slogPlot=0; else slogPlot=1; end

%% prompt user for checking imaging parameters
%reference values
kB=1.38e-23;  %Boltzman's constant (m^2 kg s^-2 K^-1)
temp=310;  %Temperature (in K)
pxlScale=0.117;  %pixel size (in μm)  (63X objective)
rad=0.1;  %radius of probe particle (in μm)
aqtime=10.57; %total time for image acquisition (must be same for all files)
nFrames=1000; %total number of frames acquired in timelapse
maxTau=800; %maximum number of frames to be averaged for MSD calculation

prompt2User = {'Temperature (in K) ?:',...
    '1 pixel = x μm, x = ?:',...
    'Particle radius (in μm) ?',...
    'Total image acquisition time (in s)'...
    'Number of frames in timelapse ?', ...
    'Maximum frames for MSD calculation (eg, 200 out of 1000 frames) ?'};
dialogBoxTitle = 'Parameters [MUST BE SAME FOR ALL FILES]';
dims = [1 100];
defaultAns = {num2str(temp),num2str(pxlScale),num2str(rad),...
    num2str(aqtime),num2str(nFrames),num2str(maxTau)};
answer = inputdlg(prompt2User,dialogBoxTitle,dims,defaultAns);

temp=str2double(cellstr(answer(1)));
pxlScale=str2double(cellstr(answer(2)));
rad=str2double(cellstr(answer(3)));
aqtime=str2double(cellstr(answer(4)));
nFrames=str2double(cellstr(answer(5)));
maxTau=str2double(cellstr(answer(6)));

%% Convert parameters to SI units
kBT=kB*temp;
pxlScale=pxlScale*(1e-6);  %pixel size (in m)
rad=rad*1e-6;  %radius of probe particle (in m)

%% begin operation
ifile=1;
while (ifile<=nf)
    infile=FileNames(ifile).name; %choose a file, name it 'infile'
    insheet=1;
    %find column with trajectory numbers, read data accordingly
    ColCheck=xlsread(infile,insheet,'A:A');
    colCheckMark=0;
    for iiCC=2:1:(length(ColCheck)-2)
        if ((ColCheck(iiCC)==ColCheck(iiCC+1))&&(ColCheck(iiCC+1)==ColCheck(iiCC+2))),
            colCheckMark=1; break,
        end
    end
    if (colCheckMark==1)
        P=xlsread(infile,insheet,'A:A');    %Particle Trajectory
        F=xlsread(infile,insheet,'B:B');    %Frame Number
        X=xlsread(infile,insheet,'C:C');    %X co-ordinate
        Y=xlsread(infile,insheet,'D:D');    %Y co-ordinate
    else
        P=xlsread(infile,insheet,'B:B');    %Particle Trajectory
        F=xlsread(infile,insheet,'C:C');    %Frame Number
        X=xlsread(infile,insheet,'D:D');    %X co-ordinate
        Y=xlsread(infile,insheet,'E:E');    %Y co-ordinate
    end
    Y=-Y;  %since Y data in image is reflected for plot image data
    
    %store trajectory numbers in TjN array
    TjN(1)=P(1); nTrj=1;
    for iP=1:1:(length(P)-1)
        if (P(iP+1)~=P(iP))
            nTrj=nTrj+1;
            TjN(nTrj)=P(iP+1);
        end
    end
    
    %store (X,Y) values for each trajectories separately
    AF=zeros(nFrames,nTrj); AX=zeros(nFrames,nTrj); AY=zeros(nFrames,nTrj);    
    for ii=1:1:length(P)
        for iTrj=1:1:nTrj
            if (TjN(iTrj)==P(ii))
                AF((F(ii)+1),iTrj)=F(ii);
                AX((F(ii)+1),iTrj)=X(ii);
                AY((F(ii)+1),iTrj)=Y(ii);
            end
        end
    end
    
    %filter out missing X,Y data and short trajectories
    %
    iTrj=1;
    while (iTrj<=nTrj)
        emptyCount=0;
        for ii=1:1:(nFrames-1),
            if ((AF(ii+1,iTrj)-AF(ii,iTrj))~=1),
                emptyCount=ii;
                break,
            end
        end
        if (emptyCount>0),
            AF(:,iTrj)=[]; AX(:,iTrj)=[]; AY(:,iTrj)=[];
            sMat=size(AF);  nTrj=sMat(2);
        end
        if (emptyCount==0), iTrj=iTrj+1; end
    end
    if (size(AF,2)==0),
        messageX = sprintf('no valid trajectory in file number %d',ifile)
        FileNames(ifile)=[]; nf=length(FileNames);
        continue,
    end
    %
        
    %to add different colours to different time points
    %
    CC=zeros(nFrames,3);
    for ic=1:1:nFrames
        CC(ic,1) = AF(ic,1)/nFrames;
        CC(ic,2) = (1-AF(ic,1)/nFrames);
        CC(ic,3) = 0;
    end
    %

    %plot raw trajectories
    %
    if (plotAllTraj==1)
    for itrj=1:1:nTrj
        figure('units','normalized','position',[.08 .25 .4 .6])
        scatter(AX(1:nFrames,itrj),-AY(1:nFrames,itrj),40,CC,'filled')
        hold on
        plot(AX(1:nFrames,itrj),-AY(1:nFrames,itrj))
        set(gca,'Color',[bR bG bB]);
        straw = sprintf('Center of Mass Trajectory');
        %straw = sprintf('Raw Trajectory of Result%d',TjN(itrj));
        title(straw)
        xlabel('\leftarrow X (\mum) \rightarrow')
        ylabel('\leftarrow Y (\mum) \rightarrow')
    end
    end
    %

    %subtracting trajectories from 'ref1' trajectory
    %
    if (subtractRefTraj==1)
    for itrj=1:1:nTrj
       
        if (itrj==ref1)
            continue
        end
    
        XX=AX(:,ref1)-AX(:,itrj);
        YY=AY(:,ref1)-AY(:,itrj);
        AX(:,itrj)=XX;
        AY(:,itrj)=YY;
    
        %plot subtracted trajectories
        %
        if (plotSubtractedTraj==1)
        figure('units','normalized','position',[.5 .25 .4 .6])
        scatter(XX,YY,5,CC,'filled')
        set(gca,'Color',[bR bG bB]);
        strsub=sprintf('Subtracted Trajectory %d-%d',TjN(itrj),ref1);
        title(strsub)
        xlabel('\leftarrow X (\mum) \rightarrow')
        ylabel('\leftarrow Y (\mum) \rightarrow')
        end
        %}
    end
    end
    %}

    %CALCULATING MSD OF TRAJECTORIES
    %
    if (maxTau>round(nFrames*0.95)), maxTau=round(nFrames*0.95); end
    PTAU=1:maxTau; %to store lag time (tau) values for plotting against MSD
    PTAU=(PTAU')*aqtime/nFrames;
    MSD=zeros(maxTau,nTrj); %to store Mean Squared Displacements
    D=zeros(maxTau,nf);  %to store Diffusion coefficients for each Tau
     
    for i=1:1:nTrj
        %not calculating MSD of reference trajectory
        if (i==ref1), continue, end
        %calculating MSD
        %
        nSum=0; %number of times summation is performed
        sumSD=0; %to store sum of squared displacements per tau value
        for tau=1:1:maxTau     %vary the value of Tau
            nSum=0; sumSD=0;
            for jj=1:1:maxTau %for MSD with sliding start points
            for j=jj:1:(nFrames-tau)
                %calculate sum of Squared Displacements for each Tau value
                sumSD=sumSD+...
                    ((AX(j+tau,i)-AX(j,i))^2)+((AY(j+tau,i)-AY(j,i))^2);
                nSum=nSum+1;
            end
            end
            %calculate the Mean Squared Displacement
            MSD(tau,i)=sumSD*(pxlScale^2)/nSum; %in m^2
            %calculate the Diffusion coefficient for each Tau
            D(tau,i)=MSD(tau,i)/(2*2*PTAU(tau)); %in m^2 s^-1
        end
    end

    %to remove the column that has all zeros (corresponding to ref1)
    MSD(:,find(sum(abs(MSD))==0))=[];
    sz=size(MSD);
    
    MSD = smoothdata(MSD,"movmean",5); %moving average smoothing
    D = smoothdata(D,"movmean",5); %moving average smoothing 
    
    lnPTAU=log(PTAU);
    lnMSD=log(MSD);
    Omega=1./PTAU; %Frequency array

    %CALCULATE G* USING THE LOCAL POWER LAW APPROXIMATION
    %
    lnPTAU2=zeros(length(lnPTAU),sz(2));%to store ln(tau) for plotting
    Alpha=zeros(length(lnPTAU),sz(2));  %to store scaling coefficients
    for i=1:1:sz(2)
        lnPTAU2(:,i)=lnPTAU;
    end
    %fit lnMSD to an nth degree polynomial
    nth=5;
    PFitCoef=zeros((nth+1),(sz(2)+1));
    gmAlpha=zeros(length(lnPTAU),sz(2));
    GlobalFit=zeros(length(lnPTAU),sz(2));
    for icr=1:1:sz(2)
        %get fitting coefficients
        PFitCoef(:,icr)=polyfit(lnPTAU2(:,icr),lnMSD(:,icr),nth);
        %generate new curve
        GlobalFit(:,icr)=polyval(PFitCoef(:,icr),lnPTAU2(:,icr));
        for k = 1:1:nth
            Alpha(:,icr)=Alpha(:,icr)+k*PFitCoef(nth-k+1,icr).*lnPTAU2(:,icr).^(k-1);
        end
        gmAlpha(:,icr)=0.457.*(1+Alpha(:,icr)).^2-1.36.*(1+Alpha(:,icr))+1.90;
        %approximation value from Mason, gmAlpha = gamma(1 + Alpha);
    end
    
    Gstar=abs(2*kBT./(3*pi*rad*MSD.*gmAlpha));
    Gelas=Gstar.*cos(pi*Alpha/2);
    Gvisc=Gstar.*sin(pi*Alpha/2);
    RelaxTime=atan(Gvisc./Gelas);
        
    %plot results
    
    MSD=MSD*1e18; D=D*1e18; %convert length from m to nm for easy plotting
    
    if (displayPlots==1)
    figure('units','normalized','position',[.0 .05 0.99 .85])
    %plot MSD vs tau in log-log scale
    subplot(2,3,1)
    loglog(PTAU,MSD,'LineWidth',2)
    ylim([5e-1 1e6]);
    if (subtractRefTraj==0)
        str1=sprintf('MSD (no drift subtraction)');
    else
        str1=sprintf('MSD (drift subtracted)');
    end
    title(str1)
    xlabel('\tau (lag time in s)')
    ylabel('MSD(nm^2)')
    %plot Diffusion rate (D) vs tau
    subplot(2,3,2)
    loglog(PTAU,D,'LineWidth',2)
    ylim([1.0 1e4]);
    title('Diffusion Rate D vs \tau');
    xlabel('\tau (lag time in s)')
    ylabel('D( nm^2/s )')
    %plot Gelas against frequency Omega
    subplot(2,3,4)
    if (logPlot==1)
        loglog(Omega(1:end),Gelas,'LineWidth',2);
        xlim([0.01 3e2]);
        ylim([5e-4 1e3]);
    end
    if (slogPlot==1)
        semilogx(Omega(1:end),Gelas,'LineWidth',2);
        xlim([0.01 3e2]);
        ylim([0 1000]);
    end
    title('G^{elastic} vs \omega (Power Law)')
    xlabel('\omega = 1/\tau (Hz)')
    ylabel('Elastic Modulus (Pa)')
    %plot Gvisc against frequency Omega
    subplot(2,3,5)
    if (logPlot==1)
        loglog(Omega(1:end),Gvisc,'LineWidth',2);
        xlim([0.01 3e2]);
        ylim([5e-4 1e3]);
    end
    if (slogPlot==1)
        semilogx(Omega(1:end),Gvisc,'LineWidth',2);
        xlim([0.01 3e2]);
        ylim([0 500]);
    end
    title('G^{viscous} vs \omega (Power Law)')
    xlabel('\omega = 1/\tau (Hz)')
    ylabel('Viscous Modulus (Pa s)')
    %plot Relaxation Time (Phase Angle) against frequency Omega
    subplot(2,3,6)
    semilogx(Omega(1:end),RelaxTime,'LineWidth',2);
    xlim([0.01 3e2]);
    ylim([-2 2]);
    title('Relaxation Time T_{\phi} vs \omega (Power Law)')
    xlabel('\omega = 1/\tau (Hz)')
    ylabel('T_{\phi} (s)')
    end

    %FIND POWER SPECTRUM OF DISPLACEMENT
    nTrj2=nTrj; %recalibrating number of trajectories to nff
    if (ref1>0)
        nTrj2=nTrj-1;
    end
    Disp=zeros(nFrames,nTrj2); FDisp=zeros(nFrames,nTrj2); Mag=zeros(nFrames,nTrj2);
    ndh=ceil(nFrames/2);    Mag1=zeros(ndh,nTrj2);
    %calculate displacements
    for i=1:1:nTrj2
        for j=2:1:nFrames
            Disp(j,i)=sqrt((AX(j,i)-AX(j-1,i))^2)+((AY(j,i)-AY(j-1,i))^2);
        end
    end
    Disp=Disp*pxlScale;
    %calculate FFT(Disp), its magnitude and phase
    for i=1:1:nTrj2
        FDisp(:,i)=fft(Disp(:,i));
        Mag(:,i)=abs(FDisp(:,i));
    end

    %calculate Single Sided Amplitude Spectrum, Mag1
    Mag=Mag/nFrames;
    Mag1(1:ndh,1:nTrj2)=Mag(1:ndh,:);  %folding upto Nyquist frequency
    Mag(:,find(sum(abs(Mag))==0))=[];
    Mag1(1,:)=0; %eliminate noise (DC component?) for better visualization
    Mag1((2:end-1),:)=2*Mag1((2:end-1),:);
    %make Frequency array
    Freq=(0:1:(ndh-1))/nFrames;
    Freq=(Freq')*nFrames/aqtime;
    %plot Single Sided Power Spectrum
    if (displayPlots==1)
    subplot(2,3,3)
    plot(Freq,Mag1(:,:),'LineWidth',1)
    ylim([0 5e-9]);
    title('PowerSpectrum-Displacement')
    xlabel('Sampling frequency (Hz)')
    ylabel('Single Sided Amplitude')
    end
    
    %store significant results
    
    if (saveMeanRheologyCurves==1)
        MSDmean(:,ifile)=mean(MSD,2);
        Dmean(:,ifile)=mean(D,2);
        GelasMean(:,ifile)=mean(Gelas,2);
        GviscMean(:,ifile)=mean(Gvisc,2);
        RelaxTimeMean(:,ifile)=mean(RelaxTime,2);
    end
    
    if (saveTwoFrequencyValues==1)
    FishOmega1Hz=abs(Omega-1); [k1min,index1Hz]=min(FishOmega1Hz);
    FishOmega10Hz=abs(Omega-10); [k10min,index10Hz]=min(FishOmega10Hz);
    for iTrj=1:1:nTrj2
        Gelas1Hz(iTrj,ifile)=Gelas(index1Hz,iTrj);
        Gelas10Hz(iTrj,ifile)=Gelas(index10Hz,iTrj);
        Gvisc1Hz(iTrj,ifile)=Gvisc(index1Hz,iTrj);
        Gvisc10Hz(iTrj,ifile)=Gvisc(index10Hz,iTrj);
        RelxT1Hz(iTrj,ifile)=RelaxTime(index1Hz,iTrj);
        RelxT10Hz(iTrj,ifile)=RelaxTime(index10Hz,iTrj);
        X1(iTrj,ifile)=AX(1,iTrj);
        Y1(iTrj,ifile)=(-1)*AY(1,iTrj);
    end
    end
    ifile=ifile+1;
end

%publish significant results
fnames=struct2cell(FileNames); fnames2=(fnames(1,:));

if (saveMeanRheologyCurves==1)
    outfileC=['MEANandSD-allFiles.xls'];
    %
    MSDmean(:,find(sum(abs(MSDmean))==0))=[];
    GelasMean(:,find(sum(abs(GelasMean))==0))=[];
    GviscMean(:,find(sum(abs(GviscMean))==0))=[];
    RelaxTimeMean(:,find(sum(abs(RelaxTimeMean))==0))=[];
    Dmean(:,find(sum(abs(Dmean))==0))=[];
    %
    l1={'Tau','mean MSD','sd MSD'};
    l2={'s','nm^2'};
    xlswrite(outfileC,l1,1,'A1');
    xlswrite(outfileC,l2,1,'A2');
    xlswrite(outfileC,PTAU,1,'A4');
    xlswrite(outfileC,mean(MSDmean,2),1,'B4');
    xlswrite(outfileC,std(MSDmean,0,2),1,'C4');
    xlswrite(outfileC,fnames2,1,'D1');
    xlswrite(outfileC,MSDmean,1,'D4');
    
    l1={'Omega','mean Gelastic','sd Gelastic'};
    l2={'Hz','Pa'};
    xlswrite(outfileC,l1,2,'A1');
    xlswrite(outfileC,l2,2,'A2');
    xlswrite(outfileC,Omega,2,'A4');
    xlswrite(outfileC,mean(GelasMean,2),2,'B4');
    xlswrite(outfileC,std(GelasMean,0,2),2,'C4');
    xlswrite(outfileC,fnames2,2,'D1');
    xlswrite(outfileC,GelasMean,2,'D4');
    
    l1={'Omega','mean Gviscous','sd Gviscous'};
    l2={'Hz','Pa s'};
    xlswrite(outfileC,l1,3,'A1');
    xlswrite(outfileC,l2,3,'A2');
    xlswrite(outfileC,Omega,3,'A4');
    xlswrite(outfileC,mean(GviscMean,2),3,'B4');
    xlswrite(outfileC,std(GviscMean,0,2),3,'C4');
    xlswrite(outfileC,fnames2,3,'D1');
    xlswrite(outfileC,GviscMean,3,'D4');
    
    l1={'Omega','mean RelaxTime','sd RelaxTime'};
    l2={'Hz','s'};
    xlswrite(outfileC,l1,4,'A1');
    xlswrite(outfileC,l2,4,'A2');
    xlswrite(outfileC,Omega,4,'A4');
    xlswrite(outfileC,mean(RelaxTimeMean,2),4,'B4');
    xlswrite(outfileC,std(RelaxTimeMean,0,2),4,'C4');
    xlswrite(outfileC,fnames2,4,'D1');
    xlswrite(outfileC,RelaxTimeMean,4,'D4');
    
    l1={'Tau','mean D','sd D'};
    l2={'s','nm^2/s'};
    xlswrite(outfileC,l1,5,'A1');
    xlswrite(outfileC,l2,5,'A2');
    xlswrite(outfileC,PTAU,5,'A4');
    xlswrite(outfileC,mean(Dmean,2),5,'B4');
    xlswrite(outfileC,std(Dmean,0,2),5,'C4');
    xlswrite(outfileC,fnames2,5,'D1');
    xlswrite(outfileC,Dmean,5,'D4');
end

if (saveTwoFrequencyValues==1)
outfileM=['Measurements1Hz10Hz.xls'];
%headings
line1st1={'Gels1Hz','Gvsc1Hz','RlxT1Hz','Gels10Hz','Gvsc10Hz','RlxT10Hz','X','Y'};
line2st1={'Pa','Pa s','s','Pa','Pa s','s'};
line1st2={'Mean','StDev','Mean','StDev','Mean','StDev','Mean','StDev',...
    'Mean','StDev','Mean','StDev'};
line2st2={'Pa','Pa','Pa s','Pa s','s','s','Pa','Pa','Pa s','Pa s','s','s'};
line3st2={'Gelas1Hz','Gelas1Hz','Gvisc1Hz','Gvisc1Hz','RelxT1Hz',...
    'RelxT1Hz','Gelas10Hz','Gelas10Hz','Gvisc10Hz','Gvisc10Hz',...
    'RelxT10Hz','RelxT10Hz'};

AvgGelas1Hz=zeros(nf,1); AvgGelas10Hz=zeros(nf,1);
StdGelas1Hz=zeros(nf,1); StdGelas10Hz=zeros(nf,1);
AvgGvisc1Hz=zeros(nf,1); AvgGvisc10Hz=zeros(nf,1);
StdGvisc1Hz=zeros(nf,1); StdGvisc10Hz=zeros(nf,1);
AvgRelxT1Hz=zeros(nf,1); AvgRelxT10Hz=zeros(nf,1);
StdRelxT1Hz=zeros(nf,1); StdRelxT10Hz=zeros(nf,1);

for ifile=1:1:nf
    %write headings
    xlswrite(outfileM,line1st1,ifile,'A1');
    xlswrite(outfileM,line2st1,ifile,'A2');
    %write values
    AllVal=cat(2,Gelas1Hz(:,ifile),Gvisc1Hz(:,ifile),RelxT1Hz(:,ifile),...
        Gelas10Hz(:,ifile),Gvisc10Hz(:,ifile),RelxT10Hz(:,ifile),X1(:,ifile),Y1(:,ifile));
    xlswrite(outfileM,AllVal,ifile,'A4');
    %remove extra 0's in the arrays
    Gelas1Hz(Gelas1Hz==0) = NaN;
    Gelas10Hz(Gelas10Hz==0) = NaN;
    Gvisc1Hz(Gvisc1Hz==0) = NaN;
    Gvisc10Hz(Gvisc10Hz==0) = NaN;
    RelxT1Hz(RelxT1Hz==0) = NaN;
    RelxT10Hz(RelxT10Hz==0) = NaN;
    %calculate Mean and Standard Deviation values
    AvgGelas1Hz(ifile)=mean(Gelas1Hz(:,ifile),'omitnan');
    AvgGelas10Hz(ifile)=mean(Gelas10Hz(:,ifile),'omitnan');
    StdGelas1Hz(ifile)=std(Gelas1Hz(:,ifile),'omitnan');
    StdGelas10Hz(ifile)=std(Gelas10Hz(:,ifile),'omitnan');    
    AvgGvisc1Hz(ifile)=mean(Gvisc1Hz(:,ifile),'omitnan');
    AvgGvisc10Hz(ifile)=mean(Gvisc10Hz(:,ifile),'omitnan');
    StdGvisc1Hz(ifile)=std(Gvisc1Hz(:,ifile),'omitnan');
    StdGvisc10Hz(ifile)=std(Gvisc10Hz(:,ifile),'omitnan');    
    AvgRelxT1Hz(ifile)=mean(RelxT1Hz(:,ifile),'omitnan');
    AvgRelxT10Hz(ifile)=mean(RelxT10Hz(:,ifile),'omitnan');
    StdRelxT1Hz(ifile)=std(RelxT1Hz(:,ifile),'omitnan');
    StdRelxT10Hz(ifile)=std(RelxT10Hz(:,ifile),'omitnan');
end
%write Mean and Standard Deviation values on the last sheet
AllStats=cat(2,AvgGelas1Hz,StdGelas1Hz,AvgGvisc1Hz,StdGvisc1Hz,...
    AvgRelxT1Hz,StdRelxT1Hz,AvgGelas10Hz,StdGelas10Hz,AvgGvisc10Hz,...
    StdGvisc10Hz,AvgRelxT10Hz,StdRelxT10Hz);
xlswrite(outfileM,line1st2,(nf+1),'A1');
xlswrite(outfileM,line2st2,(nf+1),'A2');
xlswrite(outfileM,line3st2,(nf+1),'A3');
xlswrite(outfileM,AllStats,(nf+1),'A4');
xlswrite(outfileM,fnames2',(nf+1),'N4');
end