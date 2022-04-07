%to read an image stack of FRAP data, find the nominal radius,
%and measure mean intensity within the circle of interest
%(nominal radius=radius at exp(-2) times the max intensity)
%image stack must have the filename like 't20.42s_MMStack_Pos0.ome'
%code reads '20.42s' between 't' and the first '_'
%where 20.42s is the acquisition time

clear all

%background intensity (just imaging beam without fluoresecent substances)
background=100;

FileNames=dir('t*.tif');  %array stores file names in current folder
nf=length(FileNames);   %stores number of files read

for imj=1:1:nf
    
    fname=FileNames(imj).name; %read one image stack file
    info=imfinfo(fname);%obtain all information about image stack
    nimg=numel(info);%obtain number of images in stack
    
    %read 1st image and make a binary image
    %every pixel above intensity value of exp(-2)*max intensity is white
    %all others are black
    %
    A=imread(fname,1,'Info',info);
    A=double(A);
    [nx,ny]=size(A);
    
    Abinary=zeros(nx,ny);
    %background=min(min(A))*1.0;
    
    exp2Maxima=exp(-2)*max(max(A))+(1-exp(-2))*background;
    for i=2:1:(nx-1)
        for j=2:1:(ny-1)
            if (A(i,j)>=exp2Maxima)                    
                Abinary(i,j)=1;
            end
        end
    end
    
    
    %eliminate stray bright pixels
    %
    nbr=0;
    for i=2:1:(nx-1)
        for j=2:1:(ny-1)
            if (Abinary(i,j)==1)
                nbr=0;
                for ii=(i-1):1:(i+1)
                    for jj=(j-1):1:(j+1)
               if ((Abinary(ii,jj)==1)&&((ii+jj)~=(i+j))&&((ii*jj)~=(i*j)))
                   nbr=nbr+1;
               end
                    end
                end
                if (nbr<3)
                    Abinary(i,j)=0;
                end
            end
        end
    end
    %eliminate stray dark pixels
    %
    nbr=0;
    for i=2:1:(nx-1)
        for j=2:1:(ny-1)
            if (Abinary(i,j)==0)
                nbr=0;
                for ii=(i-1):1:(i+1)
                    for jj=(j-1):1:(j+1)
               if ((Abinary(ii,jj)==0)&&((ii+jj)~=(i+j))&&((ii*jj)~=(i*j)))
                   nbr=nbr+1;
               end
                    end
                end
                if (nbr<3)
                    Abinary(i,j)=1;
                end
            end
        end
    end
    %
    
    %{
    %center of disc corresponds to point of maxima of disc
    for i=2:1:(nx-1)
        for j=2:1:(ny-1)
            if (A(i,j)==max(max(A)))                    
                centerx=i;
                centery=j;
            end
        end
    end    
    %}
    
    %find center of disc using centroid formula
    sumi=0; sumj=0; counti=0; countj=0;
    for i=2:1:(nx-1)
        for j=2:1:(ny-1)
            if (Abinary(i,j)==1) 
                sumi=sumi+i;
                sumj=sumj+j;
                counti=counti+1;
                countj=countj+1;
            end
        end
    end
    centerx=ceil(sumi/counti);
    centery=ceil(sumj/countj);
    
    %to find the nominal diameter, find the distance between the farthest
    %1s on the binary image
    %
    px1=0; px2=0;
    for j=2:1:(ny-1)
        if (Abinary(centerx,j)==1)
            if (Abinary(centerx,(j+1))~=0)
                px1=j;
                break
            end
        end
    end
    for j=(ny-1):(-1):2
        if (Abinary(centerx,j)==1)
            if (Abinary(centerx,(j-1))~=0)
                px2=j;
                break
            end
        end
    end
    radcol=(px1-px2)/2;
    py1=0; py2=0;
    for i=2:1:(nx-1)
        if (Abinary(i,centery)==1)
            if (Abinary((i+1),centery)~=0)
                py1=i;
                break
            end
        end
    end
    for i=(nx-1):(-1):2
        if (Abinary(i,centery)==1)
            if (Abinary((i-1),centery)~=0)
                py2=i;
                break                
            end
        end
    end
    
    radrow=(py2-py1)/2;
    Radn=radcol;
    if (radrow>=radcol)
        Radn=radrow;
    end
    %
    
    %paint all pixels on the nominal radius gray to check calculations
    %also, paint center grey
    %
    width=0;
    for i=2:1:(nx-1)
        for j=2:1:(ny-1)
            width=sqrt(((i-centerx)^2)+((j-centery)^2));
            if ((width>=Radn-1)&&(width<=Radn+1))                
                Abinary(i,j)=0.5;
            end
        end
    end
    Abinary(centerx,centery)=0.5;
    figure,image(Abinary*255);
    axis image
    colormap(gray(256));
    %
    
    
    %to make a video file of the image stack
    %{
    vidname=strcat(fname,'_vid.avi');
    vidObj=VideoWriter(vidname);
    vidObj.FrameRate=23;
    open(vidObj);
    %}
    
    %
    I=zeros(nimg,1);
    for iimg=1:1:nimg
        A=imread(fname,iimg,'Info',info);
        A=double(A);
        [nx,ny]=size(A);
        %(just for fun) run and save an animation of the 2.5D contour
        %{
        Img=surf(A);
        shading interp
        title('2.5D view')
        xlabel('x axis\rightarrow')
        ylabel('{\leftarrow} y axis')
        zlabel('intensity count\rightarrow')
        drawnow;
        %{
        vidframe=getframe;
        writeVideo(vidObj,vidframe);
        %}
        refreshdata(Img)
        %}
        
        %calculate intensity within imaging area inside nominal radius
        %nominal radius=radius at exp(-2) times the max intensity
        %
        intcount=0; pixcount=0;
        for i=2:1:(nx-1)
            for j=2:1:(ny-1)                    
                width=sqrt(((i-centerx)^2)+((j-centery)^2));
                if (width<Radn)
                    intcount=intcount+A(i,j);
                    pixcount=pixcount+1;
                end
            end
        end
        I(iimg)=intcount/pixcount;
        %
    end
    %close(vidObj);
    %
    %obtaining acquisition time from filenames
    aqTSEC=0; decpnt=0;
    dfname=double(fname);
    lname=length(dfname);
    for iname=1:1:lname        
        if (dfname(iname)>=48)&&(dfname(iname)<=57)   %48='0',57='9'
            aqTSEC=(10*aqTSEC)+(dfname(iname)-48);
            if (dfname(iname-decpnt-1)==46)
                decpnt=decpnt+1;
            end
            if (dfname(iname+1)==115)
                break
            end            
        end
    end
    
    aqTSEC=aqTSEC/(10^decpnt);  %total acquisition time
    nd=length(I);  %number of data points
    frametime=aqTSEC/nd;  %time period between consecutive frames
    
    %get acquisition time as string from file name
    dfname=double(fname);
    lname=length(dfname);
    for i=1:1:lname
        if ((dfname(i)>=48)&&(dfname(i)<=57))   %48='0',57='9'
            spos=i;
            break
        end
    end
    for i=1:1:lname
        if (dfname(i)==95)  %95='_'
            epos=i;
            break
        end
    end
    
    time=fname((spos-1):(epos-1));
    T=1:1:nimg; T=T';
    fileCallSign=fname((epos+1):(lname-4));
    
    %write data into a .xls file
    %
    outfilename=strcat('D:\FRAPdata\Results_',time,fileCallSign,'.xls');
    outfile0=[outfilename];
    xlswrite(outfile0,T,1,'A1');
    xlswrite(outfile0,I,1,'C1');
    %
    
    
    %initialize variables
    nd=length(I);  %number of data points
    bin=8; minbin=9999999; maxbin=0; minI=9999999; maxI=0;
    minf=0; maxf=0;
    Sumbin=zeros((nd-bin+1),1);
    
    %calculate sum of intensities per bin
    for j=1:1:(nd-bin+1)
        for k=j:1:(j+bin-1)
            Sumbin(j)=Sumbin(j)+I(k);
        end
    end
    
    %calculate prebleach intensity by averaging intensities from 3 bins    
    prebleach=(Sumbin(1)+Sumbin(2)+Sumbin(3))/(3*bin);
    
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
    stepsz=(I(blStart-1)-I(minf))/tCount;
    for j=blStart:1:(minf-1)
        tCount=tCount-1;
        I(j)=stepsz*(tCount+1)+I(minf);
    end
    %assign bleaching time
    %
    bTSEC=0.005;
    %}
    %or measure it
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
    startp=minf+300;  %curve to be fitted upto given number of data points
    %}
    lwd=3; %linewidth for all curves
    fsz=14; %fontsize for texts
    
    halfFlag=0;
    %calculating half-life of recovery
    %
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
    %
    figure('units','normalized','position',[.1 .1 .8 .8])
    if (halfFlag==1)
        plot(T,I,'b-',T,Pb,'c-',T(1:(half+minf)),Hval,'m-','linewidth',lwd)
        srt1='\tau_D';    srt2=sprintf(' = %0.3f s',tauDh(i));
        sdc=sprintf(' = %0.4f',Dhalf(i)); sdc2=' \mum^2/s,';
        srt=strcat(srt1,srt2,', D',sdc,sdc2);
    else
        plot(T,I,'b-',T,Pb,'c-','linewidth',lwd)
    end
    sat=strcat('Acquisition time = ',time);
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
    sfname = strcat(sfname0,'_Halflife.png');
    print(sfname,'-dpng')
    %}
    
end