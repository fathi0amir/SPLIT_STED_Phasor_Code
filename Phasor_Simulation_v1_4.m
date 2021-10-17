%v1_3 the curve fitting noise reduction is added.
%V1_4 added the uni-cricle underlying the 2D phasor plot
clear

%---------------------------
%Define constants
%---------------------------
iniPrm.sgnCnt=10;  %maximum counts
iniPrm.bckGrnd=1;  %Uniform backgorund level added to the image.
iniPrm.nsLvl=0.15;    %Noise level add through the poisson distribution 
iniPrm.deNs=0;     %input 1 for denoise. input 0 for noisy.
iniPrm.bmWst=167;  %Beam waist
iniPrm.tau0=2.5;   %lifetime
iniPrm.gamma0=1/iniPrm.tau0;
iniPrm.kSat=10;    %I_STED/I_Saturation Ratio
iniPrm.nCmp=2;     %number of component
%iniPrm.nHrm=nCmp-1; %number of harmonics
iniPrm.pxSz=5.2;   %Pixel Size
iniPrm.rPxSz=0.01; %Pixel size
iniPrm.imSz=64;    %Image size
iniPrm.tmMx=12.5;  %Maximum time interval
iniPrm.nTmLp=128;  %number of time laps 



%temporal steps
clear('time');
time=linspace(0,iniPrm.tmMx,iniPrm.nTmLp);

%intesity function 
intst=@(time,r1,r2, gamma0, kSat, bmWst)(exp(-gamma0*time)./(1+kSat*gamma0*time/2).*...
    (exp(-(1+kSat*gamma0*time/2)*2*r1^2/bmWst^2)-...
     exp(-(1+kSat*gamma0*time/2)*2*r2^2/bmWst^2))...
    );

%find the r when int(F, r1:r2)=int(f, 0, inf)/(number of components)
rBndTmp=0; rBnd=0; rBnd(1)=0;
for k=1:iniPrm.nCmp-1
    for i=1:floor(iniPrm.bmWst/iniPrm.rPxSz)
        r1=rBndTmp;
        r2=r1+(i-1)*iniPrm.rPxSz;
        intstInt(i)=trapz(time, intst(time,r1,r2, iniPrm.gamma0, iniPrm.kSat, iniPrm.bmWst))-...
            (1/iniPrm.nCmp)*trapz(time, intst(time,0,inf, iniPrm.gamma0, iniPrm.kSat, iniPrm.bmWst)); %#ok<*SAGROW>
    end
    [rMnVl, rMnInd]=min(abs(intstInt));
    rBnd(k+1)=rBndTmp+rMnInd*iniPrm.rPxSz;
    rBndTmp=rBnd(k+1);
end
rBnd(end+1)=inf;
clear('k', 'i', 'rBndTmp', 'r1', 'r2', 'rMnInd', 'rMnVl','intstInt')

%Generate only two particle in the center with the certain distance
%here it's 100 nm distanced
sPn=zeros(iniPrm.imSz);
sPn(iniPrm.imSz/2, iniPrm.imSz/2-floor(100/(2*iniPrm.pxSz)))=1;
sPn(iniPrm.imSz/2, iniPrm.imSz/2+floor(100/(2*iniPrm.pxSz)))=1;

% %Generate random particles. 
%Uncomment the following for random particle
% sPn=zeros(64);
% for i=1:64
%     for j=1:64
%         if rand()>0.9991
%             sPn(i,j)=1;
%         end
%     end
% end
% clear('i','j')

%count the number of particles in the simulated image and index
%their coordinates
pnInd=0;
k=1;
for i=1:64
    for j=1:64
        if sPn(i,j)==1
            pnInd(k,1)=i;
            pnInd(k,2)=j;
            k=k+1;
        end
    end
end
clear('k','i','j')

%---------------------------
%Convolution part, creates a stack of temporal spatial convolution of the
%image
%---------------------------
sPnSted=zeros(iniPrm.imSz,iniPrm.imSz,iniPrm.nTmLp);
for t=1:iniPrm.nTmLp
    for k=1:size(pnInd,1)
        for i=1:iniPrm.imSz
            for j=1:iniPrm.imSz
                r=sqrt((pnInd(k,1)-i)^2+(pnInd(k,2)-j)^2)*iniPrm.pxSz;
                gamma(i,j)=iniPrm.gamma0*(1+iniPrm.kSat*r^2/iniPrm.bmWst^2);
                sPnStedTemp(i,j,t)=...
                    iniPrm.sgnCnt*exp(-2*r^2/iniPrm.bmWst^2)*...
                    exp(-gamma(i,j)*time(t))+...
                    iniPrm.bckGrnd;
            end
        end
%         uncomment the following for the noiseless simulation
%         sPnSted(:,:,t)=sPnSted(:,:,t) + sPnStedTemp(:,:,t);
        %uncomment the following for the added noise simulation 
        sPnSted(:,:,t)=sPnSted(:,:,t) + sPnStedTemp(:,:,t)+...
            poissrnd(iniPrm.nsLvl*iniPrm.sgnCnt, [iniPrm.imSz iniPrm.imSz]);
        
    end
end
clear('sPnStedTemp','t','k','i','j','r','gamma');

%---------------------------
%fit the resulted noisy image with the the fitting function and then
%regenerate the whole image with the fitting result to avoid noise.
%---------------------------
if iniPrm.deNs
    fitFnc=@(a,time)a(1)*exp(-a(2)*time.^(a(3)))+a(4);
    lsqOptn=optimoptions('lsqcurvefit', 'display', 'none');
    lb = [];
    ub = [];
    for i=1:iniPrm.imSz
        for j=1:iniPrm.imSz
            ydata=squeeze(sPnSted(i,j,:)).';
            %       a0=[1st exp amp, 1st exp lifetime, 2nd exp amp,
            %        2nd exp lifetime, offset]
            a0=[max(ydata),0.4,1,min(ydata)];
            fitDt=lsqcurvefit(fitFnc,a0,time,ydata,lb,ub,lsqOptn);
            sPnSted(i,j,:)=reshape(fitDt(1)*exp(-fitDt(2)*time.^(fitDt(3)))+...
                fitDt(4),[1,1,size(time,2)]);
        end
        disp(' Denoising in Progress: (%)')
        disp(floor(100/iniPrm.imSz*i))
    end
    clear('a0', 'fitDt','ydata', 'i', 'j','lb','ub','lsqOptn');
end

%View the generated stack
% ViewImageStack(sPnSted)

%---------------------------
%calculate the real part and imaginary part of temeporal profile
%the 3rd dimentions indexes the higher harmonics of the fuorier transform.
%---------------------------
sPnStedRe=zeros(iniPrm.imSz,iniPrm.imSz,iniPrm.nTmLp);
sPnStedIm=zeros(iniPrm.imSz,iniPrm.imSz,iniPrm.nTmLp);
for i=1:iniPrm.imSz
    for j=1:iniPrm.imSz
        sPnStedRe(i,j,:)=real(fft(squeeze(sPnSted(i,j,:))))/...
            abs(sum(squeeze(sPnSted(i,j,:))));
        sPnStedIm(i,j,:)=imag(fft(squeeze(sPnSted(i,j,:))))/...
            abs(sum(squeeze(sPnSted(i,j,:))));
    end
end
clear('i','j','k');

%---------------------------
%Vector P made of real and imaginary parts of various harmonics 
%[g   s   gh2   sh2 ...] every element of this vector is an image in
%reciprocal space
%---------------------------
clear phsr
for i=1:floor((iniPrm.nCmp+1)/2)
    phsr(2*i-1,:,:)=sPnStedRe(:,:,i+1);
    phsr(2*i,:,:)=sPnStedIm(:,:,i+1);
end
clear('i')
if mod(iniPrm.nCmp,2)
    phsr=phsr(1:end-1, :,:);
end

%---------------------------
%building the reconstruction matrix
%Matrix M for the known number of component in the temporal dynamics
%[g1   g2   g3   ...]
%[s1   s2   s3   ...]
%[g1h2 g2h2 g3h2 ...]
%[s1h2 s2h2 s3h2 ...]
%....
%g=real part, s1=imaginary part
%---------------------------
for i=1:iniPrm.nCmp
   cmpFft(i,:)=fft(intst(time,rBnd(i),rBnd(i+1), iniPrm.gamma0, iniPrm.kSat, iniPrm.bmWst))./...
       abs(sum(intst(time,rBnd(i),rBnd(i+1), iniPrm.gamma0, iniPrm.kSat, iniPrm.bmWst)));
end
clear('i')
for i=1:floor((iniPrm.nCmp+1)/2)
    for j=1:iniPrm.nCmp
        cmpMtx(2*i-1,j)=real(cmpFft(j,i+1));
        cmpMtx(2*i,j)=imag(cmpFft(j,i+1));
    end
end
clear('i','j','cmpFft')
if mod(iniPrm.nCmp,2)
    cmpMtx=cmpMtx(1:end-1, :);
end

%---------------------------
%Temporal average 
%---------------------------
for i=1:iniPrm.imSz
    for j=1:iniPrm.imSz
        sPnStedSm(i,j)=sum(sPnSted(i,j,:));
    end
end
clear('i','j')

%---------------------------
%Reconstructing the images related to each component
%---------------------------
cmpMtxInv=inv(cmpMtx);
sPnStedRc=zeros(size(cmpMtx,1),iniPrm.imSz,iniPrm.imSz);
for k1=1:size(cmpMtx,1)
    for k2=1:size(cmpMtx,1)
        sPnStedRcTmp(k1,:,:)=cmpMtxInv(k1,k2)*squeeze(phsr(k2,:,:)).*...
            sPnStedSm(:,:);
        sPnStedRc(k1,:,:)=sPnStedRc(k1,:,:)+sPnStedRcTmp(k1,:,:);
    end
end
clear('k1','k2','sPnStedRcTemp')

%---------------------------
%Display the resconstructed image
%---------------------------
figure;
colormap(hot) % heat map
subplot(2,2,1);
imshow(squeeze(sPnStedRc(1,:,:)),[]);
title('Phasor Treatment of STED Image')
subplot(2,2,3);
plot(squeeze(sPnStedRc(1,32,:)));
set(gca, 'color',[1 1 1])

% subplot(2,2,5);
% plot(squeeze(sPnSted(32,iniPrm.imSz/2-floor(100/(2*iniPrm.pxSz)),:)))
%---------------------------
%Generating the STED Image for comparison purpose
%---------------------------
STED=sum(sPnSted,3);
subplot(2,2,2);
imshow(STED,[]);
title('STED Image')
subplot(2,2,4);
plot(squeeze(STED(32,:)));
colormap(hot) % heat map

%sPnStedRe and sPnStedIm are the real part and imaginary part of the phasor
%af the foriour transform.

%---------------------------
%Making the histogram of the phasor begins here
%---------------------------
%Determin the min and max of the phasor data points for only the first
%harmonics 
minRe=min(min(sPnStedRe(:,:,2)));
maxRe=max(max(sPnStedRe(:,:,2)));
minIm=min(min(sPnStedIm(:,:,2)));
maxIm=max(max(sPnStedIm(:,:,2)));

%reshape the data points of the first harmonics in 1D array
%and concatenate them to form m-x-2 array with first column
%as real paarts and second column as imaginary part
rePrt=reshape(sPnStedRe(:,:,2),[iniPrm.imSz^2,1]);
imPrt=reshape(sPnStedIm(:,:,2),[iniPrm.imSz^2,1]);
imRePrt=[rePrt, imPrt];

%record the intensities in the hist3 handler
intCnt = hist3(imRePrt, [64,64]);%increase the number of bin to 64by64
intCntT = intCnt';
intCntT(size(intCnt,1) + 1, size(intCnt,2) + 1) = 0; %add zero to the data points

%determine the boundaries from min to max in equal spacing
reBnd = linspace(minRe,maxRe,size(intCnt,1)+1);
imBnd = linspace(minIm,maxIm,size(intCnt,1)+1);



%the following plots the 3D histogram with the two the intensity map
%if only the "pcolor" part the executed the result is the intesity map
%only. "edgecolor', 'none'" is to remove the grids in this plot mode
%"'FaceColor','interp','CDataMode','auto'" is to make the 
%follow the intensity colormap. "'FaceAlpha',.65" is for transparency
%purpose. adding ZData to pcolor is to separate the two plot from 
%each other when put together 
figure;
hist3(imRePrt, [iniPrm.imSz,iniPrm.imSz],'FaceAlpha',.55,'edgecolor',...
    'none','FaceColor','interp','CDataMode','auto')
set(gcf,'renderer','opengl');
set(gca,  'color', [0.2 0.2 0.2]);
grid minor
set(gca,'MinorGridColor', 'g')
hold on
phrHnd=pcolor(reBnd, imBnd, intCntT);
set(phrHnd,'edgecolor','none') 
phrHnd.ZData=ones(size(intCntT)) * -max(max(intCnt));
% grid off
colormap(hot) % heat map
view(3); %to view in 3D
hold off

% figure; %This part makes the 2D phasor plot witout the uni-circle.
% phrHnd=pcolor(reBnd, imBnd, intCntT);
% set(phrHnd,'edgecolor','none') 

%Make the uni-circle and then overlay the 2D histogram of phasor plot
figure; %Make new figure
theta=0:0.001:1;
xCoor=0.5*sin(2*pi*theta);
yCoor=0.5*cos(2*pi*theta);
plot(yCoor,xCoor, '-g');%Plotting sin Vs cos
axis equal square;%Equal length fox X and Y axis
grid minor
set(gca, 'color', [0 0 0])
set(gca,'MinorGridColor', 'g',  'minorgridalpha', 1 )
hold on
intCntT(intCntT==0)=nan; %Set the zero values in the phasor as NAN
%so there won't be any zero value in the plot (for the sake of 
%transparency).
phrHnd=pcolor(reBnd, imBnd, intCntT);
set(phrHnd,'edgecolor','none') 
colormap(hot);
hold off

clear('imBnd','reBnd','imPrt','rePrt','imRePrt',...
    'intCnt','intCntT','minIm','minRe','maxRe','maxIm','phrHnd',...
    'theta', 'xCoor', 'yCoor')



