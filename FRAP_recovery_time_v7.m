clear all

%pixel calibration
pxltomicron=0.117;
%radius (0.5*FWHM) from fitting imaging beam profile to Gaussian function
Rad=16*pxltomicron;
%for Kang equation fit:
%nominal radius=half-width at exp(-2) height of imaging beam profile
Radn=15*pxltomicron;
%effective radius=radius of bleached profile
Rade=18*pxltomicron;
gm2=(Radn/Rade)^2;
%background intensity (just imaging beam without fluoresecent substances)
backgrnd=10;
outGraphTitle='EGFP in BSA 150mg per mL';
%recovery starts at? (frame number of the first frame after bleaching)
recFRAME=0; %if not known, put recFRAME=0
%number of frames (after bleaching) to be used for fitting
fitUPTO=300; %MUST BE A POSITIVE VALUE
%assign bleaching time
bTSEC=0.01;
%number of data points (leading up to bleaching frames) upto which 
%pre-bleach values will be plotted
if (recFRAME>0), iniD=recFRAME-2; else iniD=20; end

compileflag=1; %to compile and find average statistics of all the files in folder

FileNames=dir('Results*.xls');  %array stores file names in current folder
nf=length(FileNames);   %stores number of files read

for i=1:1:nf
    
    %read data from files
    infile=FileNames(i).name; %choose a file, name it 'infile'
    insheet=1;
    I=xlsread(infile,insheet,'C:C');   %mean intensity of ROI
    T=xlsread(infile,insheet,'A:A');   %frames
    
    %subtract background intensity
    I=I-backgrnd;
    
    %obtaining acquisition time from filenames
    aqTSEC=0; decpnt=0;
    name=double(infile);
    lname=length(name);
    for iname=1:1:lname        
        if (name(iname)>=48)&&(name(iname)<=57)   %48='0',57='9'
            aqTSEC=(10*aqTSEC)+(name(iname)-48);
            if (name(iname-decpnt-1)==46)
                decpnt=decpnt+1;
            end
            if (name(iname+1)==115)
                break
            end            
        end
    end
    
    aqTSEC=aqTSEC/(10^decpnt);  %total acquisition time
    nd=length(I);  %number of data points
    frametime=aqTSEC/nd;  %time period between consecutive frames
    
    %Generate acquisition time as string from filename for naming plots
    for isn=1:1:lname
        if ((name(isn)>=48)&&(name(isn)<=57))   %48='0',57='9'
            spos=isn;
            break
        end
    end
    for isn=1:1:lname
        if ((name(isn)==46)&&(name(isn+1)==120)) %46='.',120='x'
            epos=isn;
            break
        end
    end    
    sAQtime=name((spos-1):(epos-1));
    
    %initialize variables
    bin=2; minbin=9999999; maxbin=0; minI=9999999; maxI=0;
    lastf=nd; minf=0; maxf=0;
    Sumbin=zeros((nd-bin+1),1);
    
    %calculate sum of intensities per bin
    for j=1:1:(nd-bin+1)
        for k=j:1:(j+bin-1)
            Sumbin(j)=Sumbin(j)+I(k);
        end
    end
    
    %calculate prebleach intensity
    %if bleaching frame is known
    prebleach=0;
    if (recFRAME>0), for j=1:1:(recFRAME-1), prebleach=prebleach+I(j); end
    end
    prebleach=prebleach/(recFRAME-1);
    %if bleaching frame is unknown, by averaging intensities from 3 bins    
    if (recFRAME==0), prebleach=(Sumbin(1)+Sumbin(2)+Sumbin(3))/(3*bin); end
    
    %find the minimum and maximum sum, store its index in intensity array
    for j=1:1:(nd-bin)
        if (Sumbin(j)<minbin)
            minbin=Sumbin(j);
            mindex=j;
        end
        if (Sumbin(j)>maxbin)
            maxbin=Sumbin(j);
            maxdex=j;
        end
    end
    
    %find minimum intensity in the bin of minimum sum and get its index
    for j=(mindex-bin):1:(mindex+bin)
        if (I(j)<minI)
            minI=I(j);
            minf=j; %minf is the start of recovery phase
        end
    end
    if (recFRAME>0), minf=recFRAME; minI=I(minf); end
    
    %bleaching parameter
    Kbp=(1-minI/prebleach)*(1+gm2);
    %create time points in seconds    
    T=T*frametime;
    %create baseline of pre-bleach intensity
    Pb=zeros(length(T),1);
    Pb=Pb+prebleach;
    %do not show Intensity values during bleaching in graph
    tCount=0; blStart=0;
    for j=1:1:(minf-1), if (I(j)>prebleach*1.2), tCount=tCount+1; end
    end
    for j=1:1:(minf-1), if (I(j)>prebleach*1.2), blStart=j; break
        end
    end
    if (recFRAME>0), blStart=recFRAME; end
    
    stepsz=(I(blStart-1)-I(minf))/tCount;
    for j=blStart:1:(minf-1)
        tCount=tCount-1;
        I(j)=stepsz*(tCount+1)+I(minf);
    end
    
    %measure bleaching time from frames with very high intensity
    %{
    %find maximum intensity in the bin of maximum sum and get its index
    for j=(maxdex-bin):1:(maxdex+bin)
        if (I(j)>maxI)
            maxI=I(j);
            maxf=j; %maxf is the start of bleaching frame        
        end
    end
    bTSEC=frametime*(abs(maxf-minf));
    %}
    
    %calculate end of plateau before intensity decreases due to bleaching
    %
    for j=minf:1:(length(Sumbin)-bin-1)
        %find mean value of summed bins after minimum value is encountered
        meansumbin=0; nbin=0;
        for k=(j+1):1:(nd-bin-1)
            meansumbin=meansumbin+Sumbin(k);
            nbin=nbin+1;
        end
        meansumbin=meansumbin/nbin;
        %find mean of 'bin' number of Sumbin values
        refsumbin=0; nbin=0;
        for k =(j):1:(j+bin)
            refsumbin=refsumbin+Sumbin(k);
            nbin=nbin+1;
        end
        refsumbin=refsumbin/nbin;
        %check if this is greater than the mean (by some factor)
        if ((refsumbin/meansumbin)>1)
            lastf=j+bin;
            break
        end
    end
    %
    
    %calculate start of plateau (basically the point upto which the data is
    %fitted to the theoretical curve)
    %(before intensity starts to decrease due to bleaching)
    %{
    startp=minf;
    for j=minf:1:lastf
        if ((1-(I(j)/I(lastf)))<=0.0001)
            startp=j;%the closer I(j) is to I(lastf), closer to lastf we go
            break;
        end
    end
    compileflag=0;
    %}
    %or put any value of 'startp' as required
    %
    startp=minf+fitUPTO;  %curve to be fitted with fitUPTO frames
    %}
    lwd=3; %linewidth for all curves
    fsz=14; %fontsize for texts
    
    halfFlag=0;
    %calculating half-life of recovery
    %{
    %initialize variables
    halfFlag=1;
    diff=0; mindiff=9999999; half=0;
    %define lower & upper cutoff, and 1/2 of end point of plateau intensity
    %to get 1/2life when intensity becomes 1/2 of plateau while recovering
    lcut=(I(lastf)+minI)/2.5;
    cntr=(I(lastf)+minI)/2.0;
    ucut=(I(lastf)+minI)/1.5;
    
    for j=minf:1:lastf
        %calculate halflife time when intensity becomes half of plateau
        if (I(j)>lcut)&&(I(j)<ucut)
            diff=abs(I(j)-cntr);
            %minimize difference between 1/2 of last data point with
            %actual 1/2 value of plateau to get the best half life value 
            if mindiff>diff
                mindiff=diff;
                half=j;
            end
        end
    end
    
    Hval=zeros((half+minf),1);
    Hval=Hval+I(half); %to create an indicator line for half intensity
    halflife=frametime*(half-minf);  %half-life time
    tauDh(i)=halflife/(log(2));  %characteristic recovery time
    %Dhalf(i)=0.88*(Rad^2)/(4*halflife); %find Axelrod diffusion constant
    %Dhalf(i)=(Radn^2)/(4*tauDh(i));
    Dhalf(i)=(Radn^2)/(4*tauDh(i)-(4*bTSEC/pi));
    %}
    
    %plot intensity vs time curve with/without half-life marker
    %{
    figure('units','normalized','position',[.1 .1 .8 .8])
    if (halfFlag==1)
        plot(T,I,'b-',T,Pb,'c-',T(1:(half+minf)),Hval,'m-','linewidth',lwd)
        srt1='\tau_D';    srt2=sprintf(' = %0.3f s',tauDh(i));
        sdc=sprintf(' = %0.4f',Dhalf(i)); sdc2=' \mum^2/s,';
        srt=strcat(srt1,srt2,', D',sdc,sdc2);
    else
        plot(T,I,'b-',T,Pb,'c-','linewidth',lwd)
    end
    sat=sprintf('(Acquisition time = %0.2f s)',aqTSEC);
    stl=strcat('Intensity vs Time plot',' ',sat);
    title(stl);
    xlabel('time (in seconds)');
    ylabel('mean intensity (grayvalue) per unit area');
    dimint=[0.4 0.101 0.3 0.3];
    sbt=sprintf(' Bleaching time = %0.4f s,',bTSEC);
    sdn1='blue = total data, cyan = pre-bleach intensity baseline';
    sdn2='magenta = half-intensity indicator';
    if (halfFlag==1), strb={sbt,' ',srt,sdn1,sdn2};
    else, strb={sbt,' ',sdn1};
    end
    annotation('textbox',dimint,'String',strb,'FitBoxToText','on');
    %}
    %{
    sfname0=sprintf('t%0.2fs',aqTSEC);
    sfname = strcat(sfname0,'-Halflife.png');
    print(sfname,'-dpng')
    %}
    
    %write all the data in .xls file
    %{
    outfile0=['DataFRAP.xls'];    
    xlswrite(outfile0,T,i,'A1');
    xlswrite(outfile0,I,i,'B1');
    xlswrite(outfile0,aqTSEC,i,'C1');
    xlswrite(outfile0,bTSEC,i,'C2');
    xlswrite(outfile0,Dhalf(i),i,'C3');
    %}
    
    %selecting data during recovery time only for recovery curve fitting
    %
    Ir=zeros((startp-minf+1),1);
    Tr=zeros((startp-minf+1),1);
    Tr2=zeros((startp-minf+1),1);
    k=1;
    for j=(minf):1:(startp)
        Ir(k)=I(j);
        Tr2(k)=T(j)-T(minf);
        Tr(k)=T(j);
        k=k+1;
    end
    lastr=k-1;
    
    %points for plotting
    T2=T(minf-iniD:startp)-T(minf); I2=(I(minf-iniD:startp))/prebleach;
    
    %normalize the recovery data only
    Ir2=Ir/prebleach;
    %Ir2=Ir/Ir(startp-minf);
    %Ir2=Ir2-Ir2(1);
    %}
    
    %write recovery data in .xls file
    %{
    outfile1=['RecoveryData.xls'];
    xlswrite(outfile1,Tr2,i,'A1');
    xlswrite(outfile1,Ir2,i,'B1');
    xlswrite(outfile1,aqTSEC,i,'C1');
    xlswrite(outfile1,bTSEC,i,'C2');
    xlswrite(outfile1,Dhalf(i),i,'C3');
    %}
    
    %plot intensity vs time curve
    %{
    figure('units','normalized','position',[.04 .25 .45 .6])
    plot(T,I,'g-',Tr,Ir,'b-',T,Pb,'c-','linewidth',lwd)
    sat=sprintf('(Acquisition time = %0.2f s)',aqTSEC);
    stl=strcat('Intensity vs Time plot',' ',sat);
    title(stl);
    xlabel('time (in seconds)','fontweight','bold');
    ylabel('mean intensity','fontweight','bold');
    dimint=[0.25 0.05 0.3 0.3];    
    sbt=sprintf(' Bleaching time = %0.4f s,',bTSEC);
    sdn1='blue = recovery data, green = total data';
    sdn2='cyan = pre-bleach intensity baseline';
    strb={sbt,sdn1,sdn2};
    annotation('textbox',dimint,'String',strb,'FitBoxToText','on');
    %}
    %{
    sfname0=sprintf('t%0.2fs',aqTSEC);
    sfname = strcat(sfname0,'-IvsT.png');
    print(sfname,'-dpng')
    %}
    
    %fitting data to Kang's simplified equation
    %
    %{
    the simplified Kang equation:
    F(t) = Fi*Mf{1-Kbp/(1+gm2+AA*t)}+Fo*(1-Mf)
    where AA=2/tauD, gm2=(Radn/Rade)^2,Mf=mobile fraction
    Fi=prebleach intensity(=1 for normalized data),
    Fo=intensity just after bleaching.
    %}
    
    %define fitting function
    %fitfunc=@(p,x) (p(2)*(1-Kbp./(1+gm2+8*p(1)*Tr2/Rade^2))+(1-p(2))*Ir2(1)).*exp(-0.0065.*Tr2);
    fitfunc=@(p,x) (p(2)*(1-Kbp./(1+gm2+8*p(1)*Tr2/Rade^2))+(1-p(2))*Ir2(1));
    p=nlinfit(Tr2,Ir2,fitfunc,[1 1]);
    %initial guess  D=1 and Mf=1
    
    %The fitting parameters are D and Mf
    D(i)=p(1);   Mf(i)=p(2);
    tauD=(Rade^2)/(4*D(i)); %recovery time
    %get the best fitting curve
    %fke=(Mf(i)*(1-Kbp./(1+gm2+8*D(i)*Tr2/Rade^2))+(1-Mf(i))*Ir2(1)).*exp(-0.0065.*Tr2);
    fke=(Mf(i)*(1-Kbp./(1+gm2+8*D(i)*Tr2/Rade^2))+(1-Mf(i))*Ir2(1));
    
    FittedData=fke;
    stmf=sprintf('Mobile Fraction = %0.2f',Mf(i));
    sft='red- F(t)=F_i*M_f*(1-K_b_p/(1+\gamma^2+8*D*t/r_e^2))+F_o*(1-M_f)';
    sfname = strcat(sAQtime,'-KangEquationFit');
    %}
    
    
    doubleRateFlag=0;
    %fitting data to (A+B*exp(-C*x)+D*exp(-E*x))
    %{
    %A = Mobile Fraction; C,E = recovery rates;
    doubleRateFlag=1;
    %fit the 'double exponent inversion recovery' curve
    ft=fittype('a+b*exp(-c*x)+d*exp(-e*x)');
    spta=Ir2(lastr)-Ir2(1); sptb=Ir2(1)/2; sptc=1; sptd=Ir2(1)/2; spte=1;
    abcde=[spta;sptb;sptc;sptd;spte];
    fe2=fit(Tr2(:),Ir2(:),ft,'StartPoint',abcde);
    cValues(:,i)=coeffvalues(fe2); %store a,b,c in vector cValues
    %calculate characteristic diffusion times
    tauD1=1/abs(cValues(3,i));
    tauD2=1/abs(cValues(5,i));
    %calculate Mobile fraction
    Mf(i)=cValues(1,i);
    %store other coefficients
    sB=cValues(2,i); sD=cValues(4,i);
    
    %calculating diffusion constant 1
    %D1(i)=(Rad^2)/(4*tauD1);
    D1(i)=(Rad^2)/(4*tauD1-(4*bTSEC/pi));
    %calculating diffusion constant 2
    %D2(i)=(Rad^2)/(4*tauD2);
    D2(i)=(Rad^2)/(4*tauD2-(4*bTSEC/pi));
    
    %selecting faster components (for later use)
    if (tauD1<tauD2)
        tauD=tauD1;
        D(i)=D1(i);
        slowTau(i)=tauD2;
        slowD(i)=D2(i);
    else
        tauD=tauD2;
        D(i)=D2(i);
        slowTau(i)=tauD1;
        slowD(i)=D1(i);
    end
    
    FittedData=fe2(Tr2);
    stmf=sprintf('Mobile Fraction=%0.2f',Mf(i));
    sft='red- y=A+B*exp(-C*x)+D*exp(-E*x)';
    sfname0=sprintf('t%0.2fs',aqTSEC);
    sfname = strcat(sAQtime,'-2expFit');
    %}
    
    %OTHER FITTING MODELS
    %{
    
    
    %fitting data to y=A+B*(1-exp(-C*(t-D)))
    %{
    %{
    in the fitting equation:
    A=(Mobile Fraction-Background), B=Intensity after bleaching,
    C=Recovery Rate, D=X-axis offset.
    %}
    %define fitting function
    fitfunc=@(p,x) (p(1)+p(2)*(1-exp(-p(3)*(Tr2-p(4)))));
    sptA=(Ir2(lastr)-Ir2(1)); sptB=Ir2(1);  sptC=1;  sptD=0;
    p=nlinfit(Tr2,Ir2,fitfunc,[sptA sptB sptC sptD]);
    %initial guess
        
    %get the best fitting curve
    bfe1=(p(1)+p(2)*(1-exp(-p(3)*(Tr2-p(4)))));
    
    %calculate recovery time and diffusion coefficient
    tauD=1/p(3); Mf(i)=p(1)+p(2);
    D(i)=(Rad^2)/(4*tauD-(4*bTSEC/pi));
    
    FittedData=bfe1;
    stmf=sprintf('Mobile Fraction = %0.2f',Mf(i));
    sft='red- y=A+B*(1-exp(-C*(t-D)))';
    sfname = strcat(sAQtime,'-1expFit');
    %}
    
    
    %fitting data to Ellenberg (1997)
    %{
    %{
    y=A*(1-1/sqrt(1+(4*pi*B/Rade^2)*x))
    in the fitting equation:
    A=Mobile Fraction, B=Recovery Rate.
    %}
    %define fitting function
    fitfunc=@(p,x) p(1)*(1-1./sqrt(1+(4*pi*p(2)/Rade^2)*Tr2));
    p=nlinfit(Tr2,Ir2,fitfunc,[1 0.1]);
    %initial guess  A=1, B=0.1
        
    %get the best fitting curve
    bem=p(1)*(1-1./sqrt(1+(4*pi*p(2)/Rade^2)*Tr2));
    
    %calculate recovery time and diffusion coefficient
    tauD=1/p(2);
    D(i)=(Rad^2)/(4*tauD-(4*bTSEC/pi));
    Mf(i)=p(1);
    
    FittedData=bem;
    stmf=sprintf('Mobile Fraction = %0.2f',Mf(i));
    sft='red- y=A*(1-1/(x*sqrt(1+(4*\pi*B/r^2)))';
    sfname = strcat(sAQtime,'-Ellenberg1997Fit');
    %}
    
    %fitting data to Soumpassis Equation
    %{
    ft=fittype('a+b*exp(-c/x)*(besseli(0,c/x)+besseli(1,c/x))');
    spta=Ir2(1); sptb=(Ir2(lastr)-Ir2(1)); sptc=1;
    abc=[spta;sptb;sptc];
    fse=fit(Tr,Ir2,ft,'StartPoint',abc);
    cfValues(:,i)=coeffvalues(fse); %store a,b,c in matrix cValues
    tauD=1/(cfValues(3,i));  %calculate characteristic diffusion time
    %calculating diffusion constant
    %D(i)=(Rad^2)/(4*tauD);
    D(i)=(Rad^2)/(4*tauD-(4*bTSEC/pi));
    Mf(i)=(cfValues(1,i)+cfValues(2,i));
    
    FittedData=fse(Tr);
    stmf=sprintf('Mobile Fraction = %0.2f',Mf(i));
    sft='red- y=a+b*exp(-c/x)*(besseli(0,c/x)+besseli(1,c/x))';
    sfname = strcat(sAQtime,'-SoumpassisFit');
    %}
    
    
    %fitting data to Kang's series summation (not working)
    %{
    %{
    the Kang series summation:    
    f(t) = SUM <m=0 to inf> (-Kbp)^m/(m![1+m*(AA*t+gm2)])
    where AA=2/tauD, gm2=(Radn/Rade)^2
    then do: r(j) = Ir2(j)- f(t)
    then do: d[r(j)]/d[AA] to get the Jacobian matrix- Jacb
    for Gauss-Newton fitting
    %}    
    Jacb=zeros(lastr,1);
    Rdif=zeros(lastr,1);
    Kfit=zeros(lastr,1);
    AA=2;  %initial guess of recovery time
    runsz=100; mlim=10;
    Kbp=(Ir2(1));
    for ii=1:1:runsz   
        for j=1:1:lastr            
        jacterm=-Kbp*AA*Tr2(1)/((1+gm2+AA*Tr2(1))^2);
        kterm=1;
            for m=1:1:mlim
                mfact=1;
                for mm=1:1:m
                    mfact=mfact*mm;
                end
                Knum=((-Kbp)^m)/mfact;
                Kden=1+m*gm2+m*AA*Tr2(j);
                kterm=kterm+(Knum/Kden);
                Rnum=(Tr2(j)*AA*((-Kbp)^(m+1)))/(mfact*(m+1));
                Rden=(1+(m+1)*(AA*Tr2(j)+gm2))^2;
                jacterm=jacterm+(Rnum/Rden);
            end
            Jacb(j)=jacterm;
            Rdif(j)=Ir2(j)-kterm;
            Kfit(j)=kterm;
        end
        JacbT=Jacb';
        JtJ=JacbT*Jacb;
        AA=AA+(inv(JtJ))*JacbT*Rdif;
    end
    
    tauD=2/AA;
    %calculating diffusion constant
    %D(i)=(Rad^2)/(4*tauD);
    D(i)=(Rad^2)/(4*tauD-(4*bTSEC/pi));
    Mf(i)=max(Kfit);
    
    FittedData=Kfit;
    stmf=sprintf('Mobile Fraction = %0.2f',Mf(i));
    sft='red- f(t) = SUM <m=0 to inf> (-Kbp)^m/(m![1+m*(AA*t+gm2)])';
    sfname = strcat(sAQtime,'-KangSeriesFit');
    %}
    
    %}
    
    %separately generate Bleaching Time vs Recovery Time data
    %{
    BTs(i)=bTSEC;
    RTs(i)=tauD;
    %}
    
    
    %calculate R-squared value and plot recovery curve with fit
    %{
    R2=0; ydatamean=0; numofdata=length(Ir2); SSE=0; yvar=0;
    for j=1:1:numofdata     %calculate mean of y-data
        ydatamean=ydatamean+Ir2(j);
    end
    ydatamean=ydatamean/numofdata;
    ERR=zeros(numofdata,1);
    for j=1:1:numofdata
        ERR(j)=Ir2(j)-FittedData(j); %y(i)-f(y(i))
        SSE=SSE+ERR(j)^2; %sum of (y(i)-f(y(i)))^2
        yvar=yvar+(Ir2(j)-ydatamean)^2;   %sum of (y(i)-y-mean)^2        
    end
    R2=1-(SSE/yvar);
    
    %plot the fitted recovery curve
    figure('units','normalized','position',[.5 .2 .45 .7])
    axes('Position', [0.1 0.25 0.8 0.65]);
    plot(T2,I2,'b-',Tr2,FittedData,'r-','linewidth',lwd)
    set(gca,'fontweight','bold','fontsize',fsz);
    set(gca,'XLim',[T2(1) 1.05*T2(length(T2))])
    set(gca,'YLim',[0.0 1.2])
    set(gca,'YTick',(0.0:0.2:1.2))
    title(sfname,'fontweight','bold','fontsize',fsz);
    %xlabel('time (in seconds)','fontweight','bold','fontsize',fsz);
    set(gca,'LineWidth',2,'XTickLabel',[]);
    ylabel('Normalised mean intensity','fontweight','bold','fontsize',fsz);
    dimfit=[0.3 0.42 0.1 0.1];
    sbt= sprintf('Bleached for = %0.3f s,',bTSEC);
    sat=sprintf('(Data acquired for = %0.2f s)',aqTSEC);
    sbat=strcat(sbt,sat);
    sR2=sprintf('SSE = %0.3f, R^2 = %0.3f',SSE,R2);
    if (doubleRateFlag==0)
        srt11='\tau_D';    srt12=sprintf(' = %0.2f s',tauD);
        sdc1=sprintf(', D = %0.2f ',D(i));    sdc2 = ' \mum^2/s';
        srt=strcat(srt11,srt12,sdc1,sdc2);
        strf={srt,stmf,'blue- data',sft,sbat,sR2};
    end
    if (doubleRateFlag==1)
        srt11='\tau_{D1}';    srt12=sprintf(' = %0.3f s',tauD1);
        srt21='\tau_{D2}';    srt22=sprintf(' = %0.3f s',tauD2);
        sdc1=sprintf(', D1 = %0.3f ',D1(i));
        sdc2=sprintf(', D2 = %0.3f ',D2(i));
        sdcu='\mum^2/s, ';
        srtt1=strcat(srt11,srt12,sdc1,sdcu);
        srtt2=strcat(srt21,srt22,sdc2,sdcu);
        strf={srtt1,srtt2,stmf,'blue- data',sft,sbat,sR2};
    end
    anno=annotation('textbox',dimfit,'String',strf,'FitBoxToText','on');
    set(anno,'units','normalized','fontweight','bold','fontsize',11);
    set(anno,'EdgeColor','none','VerticalAlignment','middle');
    
    %plot ERR vs time points
    axes('Position', [0.1 0.1 0.8 0.1]);
    zerolineERR=zeros(1,numofdata);
    plot(Tr2,ERR,'b',Tr2,zerolineERR,'r','linewidth',lwd);
    set(gca,'LineWidth',2,'fontweight','bold','fontsize',fsz);
    xlabel('time (in seconds)','fontweight','bold','fontsize',fsz);
        
    %sfname1 = strcat(sfname,'.bmp');
    %print(sfname1,'-dbmp')    
    %}
    
    rTSEC(i)=tauD; %store recovery time from each file
    
    %compile recovery data to put in one spreadsheet
    %
    if (compileflag==1)
        RecoveryTime=Tr2;
        RecoveryTim2=T2;
        RecoveryData(:,i)=I2;
        RecoveryFitD(:,i)=FittedData;
    end
    %
end

%sort the bleaching time data in ascending order and compare to recovery
%times to see how bleaching time affects recovery time
%{
minbt=9999; minbtpos=0; tempbtrt=0;
for i=1:1:(nf-1)
    minbt=9999;
    for j=i:1:nf
        if (BTs(j)<=minbt)
            minbt=BTs(j);
            minbtpos=j;
        end
    end
    tempbtrt=BTs(minbtpos);
    BTs(minbtpos)=BTs(i);
    BTs(i)=tempbtrt;
    tempbtrt=RTs(minbtpos);
    RTs(minbtpos)=RTs(i);
    RTs(i)=tempbtrt;
end
%
BTs=BTs';  RTs=RTs';
%print generated data in file
%
outfile2=['bleachVSrecovery.xls'];    
xlswrite(outfile2,BTs,1,'A1');
xlswrite(outfile2,RTs,1,'B1');
%}

%Calculate Mean and Standard Deviation of recovery time (from halftime)
%{
muT=mean(tauDh);
sdT=sqrt(var(tauDh));
avgRTh=strcat(num2str(muT),'±',num2str(sdT))
meanDhalf=mean(Dhalf);
sdDhalf=sqrt(var(Dhalf));
avgDhalf=strcat(num2str(meanDhalf),'±',num2str(sdDhalf))
%}

%Calculate Mean and Standard Deviation of recovery time (from fitting)
%
muT=mean(rTSEC);
sdT=sqrt(var(rTSEC));
avgRT=strcat(sprintf('%0.2f',muT),'±',sprintf('%0.2f',sdT),' s')
meanD=mean(D);
sdD=sqrt(var(D));
avgD=strcat(sprintf('%0.2f',meanD),'±',sprintf('%0.2f',sdD),' µm^2/s')

if (doubleRateFlag==1)
    muT2=mean(slowTau);
    sdT2=sqrt(var(slowTau));
    avRT2=strcat(sprintf('%0.2f',muT2),'±',sprintf('%0.2f',sdT2),' s')
    meanD2=mean(slowD);
    sdD2=sqrt(var(slowD));
    avD2=strcat(sprintf('%0.2f',meanD2),'±',sprintf('%0.2f',sdD2),' µm^2/s')
end
avgMf=mean(Mf);

%write D and Mobile Fraction values in a spreadsheet
outfile3=strcat('DandMF-',outGraphTitle,'.xls');
D=D'; Mf=Mf'; fnames=struct2cell(FileNames); fnames2=(fnames(1,:))';
strD={'D µm^2/s'}; strMf={'Mobile Fraction'};
xlswrite(outfile3,strD,1,'A1');
xlswrite(outfile3,D,1,'A2');
xlswrite(outfile3,strMf,1,'C1');
xlswrite(outfile3,Mf,1,'C2');
if (doubleRateFlag==1)
    slowD=slowD'; strSlowD={'D_2 µm^2/s'};
    xlswrite(outfile3,strSlowD,1,'B1');
    xlswrite(outfile3,slowD,1,'B2');
end
xlswrite(outfile3,fnames2,1,'E2');
%}

%print recovery data with fitting data and the averaged data in excel files
%column 1 for time
%columns with even numbers for intensities
%columns with odd numbers>1 for fitted data
%
if (compileflag==1)
    outfile4=strcat('RecoveryDataAll-',outGraphTitle,'.xls');
    startCell=strcat('C',num2str(iniD+1));
    for isheet=1:1:nf
        xlswrite(outfile4,RecoveryTim2,isheet,'A1');
        xlswrite(outfile4,RecoveryData(:,isheet),isheet,'B1');
        xlswrite(outfile4,RecoveryFitD(:,isheet),isheet,startCell);
    end
    %print average recovery data in another spreadsheet
    [m8,n8]=size(RecoveryData);
    AvgRecData=zeros(m8,1);
    for jjj=1:1:m8, AvgRecData(jjj)=mean(RecoveryData(jjj,:)); end
    [m9,n9]=size(RecoveryFitD);
    AvgFitData=zeros(m9,1); StdFitData=zeros(m9,1);
    for jjj=1:1:m9, AvgFitData(jjj)=mean(RecoveryFitD(jjj,:)); end
    for jjj=1:1:m9, StdFitData(jjj)=sqrt(var(RecoveryFitD(jjj,:))); end
    outfile5=strcat('AverageRecoveryData-',outGraphTitle,'.xls');
    xlswrite(outfile5,RecoveryTim2,1,'A1');
    xlswrite(outfile5,AvgRecData,1,'B1');
    xlswrite(outfile5,AvgFitData,1,startCell);
    startCell2=strcat('D',num2str(iniD+1));
    xlswrite(outfile5,StdFitData,1,startCell2);
end
%

%plot recovery data in one place
%
sRA=size(RecoveryFitD);
if (compileflag==1)
    Iav=zeros(length(Ir2),1);
    Tav=RecoveryTime;    Tar=RecoveryTim2;
    for icol=1:1:nf
        Iav=Iav+RecoveryFitD(:,icol);
    end
    Iav=Iav/nf;
    figure('units','normalized','position',[.5 .2 .45 .6])
    for icol=1:1:nf
        plot(Tar,RecoveryData(:,icol),'color',[0.7,0.7,0.7],'linewidth',lwd)
        ylim([0.0 1.0])
        hold on
    end
    for icol=1:1:nf
        %plot(Tav,RecoveryAll(:,icol),'color',[0.0,0.0,0.0],'linewidth',2)
        plot(Tav,RecoveryFitD(:,icol),'r-','linewidth',2)
        ylim([0.0 1.0])
        hold on
    end
    for icol=1:1:nf
        txtn=sprintf('%d',icol);
        text(Tav((sRA(1)-1),1),RecoveryFitD((sRA(1)-1),icol),txtn)
        hold on
    end
    hold off
    %plot details
    stmf=sprintf('Average Mobile Fraction = %0.2f',avgMf);
    set(gca,'fontweight','bold','fontsize',fsz);
    set(gca,'XLim',[RecoveryTim2(1) 1.05*RecoveryTim2(length(RecoveryTim2))])
    set(gca,'YLim',[0.0 1.2])
    set(gca,'YTick',(0.0:0.2:1.2))
    %title('Recovery Data with Fit','fontweight','bold','fontsize',fsz); ctitle=0;
    customTitle=outGraphTitle; ctitle=1;
    title(customTitle,'fontweight','bold','fontsize',fsz);
    xlabel('time (in seconds)','fontweight','bold','fontsize',fsz);
    ylabel('Normalized mean intensity','fontweight','bold','fontsize',fsz);
    dimfit=[0.3 0.42 0.1 0.1];
    sbt= sprintf('Bleached for = %0.3f s',bTSEC);
    sat=sprintf(', Data acquired for = %0.2f s',aqTSEC);
    sbat=strcat(sbt,sat);
    if (doubleRateFlag==0)
        srt1=strcat('\tau_D = ',avgRT);
        sdc1=strcat(', D = ',avgD);
        srt=strcat(srt1,sdc1);
        strf={srt,stmf,'gray- data',sft,sbat};
    end
    if (doubleRateFlag==1)
        srt1=strcat('\tau_D_1 = ',avgRT); srt2=strcat('\tau_D_2 = ',avRT2);
        sdc1=strcat(', D_1 = ',avgD);     sdc2=strcat(', D_2 = ',avD2);
        srtt1=strcat(srt1,sdc1); srtt2=strcat(srt2,sdc2);
        strf={srtt1,srtt2,stmf,'gray- data',sft,sbat};
    end
    anno=annotation('textbox',dimfit,'String',strf,'FitBoxToText','on');
    set(anno,'units','normalized','fontweight','bold','fontsize',11);
    set(anno,'EdgeColor','none','VerticalAlignment','middle');
    
    sfname1 = strcat(sfname,'-ALL.bmp');
    if (ctitle==1), sfname1 = strcat(customTitle,'.bmp'); end
    print(sfname1,'-dbmp')
end
%}