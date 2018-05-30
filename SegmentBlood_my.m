function outputParams=SegmentBlood(inputParams)
% Function that reads a multi-frame volume of perfusion data
% and returns labelmap of voxels with high peak (blood)
%
% Parameters:
% inputParams.input4Dvolume: input image filename
% inputParams.inputEXROImap (optional)
% inputParams.inputINCROImap
% inputParams.outputlabelmap (output label structure, result of the processing)
% inputParams.tesla
% inputParams.r1
% inputParams.bloodlabelvalue
% outputParams.numpixelsROI
% outputParams.numpixelsAIF
% outputParams.FT
% outputParams.numframes
% outputParams.FA
% outputParams.TR
%------------------------------------------------
% Area of blood pixels to choose for AIF label
switch inputParams.animal
case 'Human'
areaAIF=50; % in mm^2; determines number of pixels used
case 'Rat'
areaAIF=2;
case 'Mouse'
%areaAIF=0.5; %For 20 piexels
areaAIF=0.5/1; %J% Only from center slice (z=3) 
end
perfvol=cli_imageread(inputParams.input4Dvolume);
perfdata=double(perfvol.pixelData);
sizeperf=size(perfdata);
eval(['t=[',perfvol.metaData.MultiVolume_FrameLabels,'];']); % time points of acquisition
t=t'/1000; % seconds
outputParams.FT=t(2)-t(1);
numframes=length(t);
outputParams.numframes=numframes;
% one of the two ROI's used as Blood label template, so at least one must
% be given
if isfield(inputParams,'inputEXROImap')
    labelvol=cli_imageread(inputParams.inputEXROImap);
    exROI=find(~labelvol.pixelData(:)); % indices of voxels to include
    if isfield(inputParams,'inputINCROImap')
        labelvol=cli_imageread(inputParams.inputINCROImap);
        incROI=intersect(exROI,find(labelvol.pixelData(:))); % indices of voxels to include
    else
        incROI=exROI;
    end
else
    if isfield(inputParams,'inputINCROImap')
        labelvol=cli_imageread(inputParams.inputINCROImap);
        incROI=find(labelvol.pixelData(:)); % area to include
    else
        disp('At least one ROI is required')
    end
end
outputParams.numpixelsROI=length(incROI);
eval(['FA=',perfvol.metaData.MultiVolume_DICOM_FlipAngle,';']);
eval(['TR=',perfvol.metaData.MultiVolume_DICOM_RepetitionTime,';']);
outputParams.FA=FA;
outputParams.TR=TR;
switch inputParams.tesla
case 1.5
T1blood=1.440; %sec
case 3
T1blood=1.630; %sec
case 7
T1blood=1.8;
end
r1=inputParams.r1;

qt=round(size(perfdata,1)/4); %J% size(perfdata,1) == 90 (Num of frames);
sigbase=squeeze(mean(perfdata(1:4,incROI))); %J% assume that frame 1-4 is before BOT
maxsig=squeeze(max(perfdata(1:qt,incROI))); % max signal in 1st quarter
meansig=squeeze(mean(perfdata(qt:2*qt,incROI))); % in 2nd quarter
stdi=squeeze(std(perfdata(2*qt:end,incROI))); % variation in 2nd half
%firstpts=squeeze(mean(perfdata(1:3,:,:,:)));
%mini=squeeze(min(perfdata));
%diffi=squeeze(mean(diff(perfdata)));
%absdiffi=squeeze(mean(abs(diff(perfdata))));
res=nrrd2res(perfvol.metaData.space_directions); % spatial resolution [mm]
bloodpixnum=round(areaAIF/(res(1)*res(2)));
%outputParams.numpixelsAIF=bloodpixnum;
outputParams.numpixelsAIF=bloodpixnum; 

bloodCharacteristics=(maxsig-sigbase).*(maxsig-meansig)./nonzero(stdi);

peaksharpness=sort(bloodCharacteristics); % NaNs at end
[~,I]=max(peaksharpness); %J% I==number of pixels in IncROI

%**********Jihun added*************%
j = 0;
for i=1:I
    if maxsig(1,i)<sigbase(1,i)
        if maxsig(1,i)<meansig(1,i)
            j = j+1;
        end
    end
end
fprintf('Number of miss pixel is %d', j);

if j>0
    %Returns error if minus*minus pixel included in bloodCharacteristics.
    disp('Psudo error flag %d', j); %just stop running.
end
%***********************************%

arteryindex=peaksharpness(I-bloodpixnum); %J% peaksharpness value of (I-20)
bloodindex = find(bloodCharacteristics > arteryindex); %J% Find index larger than arteryindex (top 20th)

%bloodi=zeros(sizeperf(2:end)); %J% Old code
bloodi=zeros(size(labelvol.pixelData)); %J% New code

bloodi(incROI(bloodindex))=1;
labelvol.pixelData = inputParams.bloodlabelvalue * bloodi;

%labelvol.pixelData = inputParams.bloodlabelvalue * bloodi;
%-------- write the .nrrd file of the blood voxels ROI
cli_imagewrite(inputParams.outputlabelmap,labelvol);

% get rid of nans, imaginary numbers, and negative numbers
%perfconc=perfconc.*~isnan(perfconc).*(perfconc==real(perfconc)).*(perfconc>0);
% assumes max blood signal is in first quarter

% Plot and Fit the AIF to a sigmoid*biexponential
%bloodi4D=permute(repmat(bloodi,[1,1,1,numframes]),[4,1,2,3]);
%Cp=squeeze(sum(sum(sum(bloodi4D.*perfconc,4),3),2))/bloodvoxnum;
hct=0.45;
bloodsignal=mean(perfdata(:,incROI(bloodindex)),2);
Cp=(1/(1-hct))*signal2conc(bloodsignal,T1blood,TR/1000,FA,r1);
if isfield(inputParams,'filename')
    savefilename=inputParams.filename;
else
    savefilename='SegmentBloodParams';
end
AIFparams=fitAIF(Cp,t,savefilename); % also plots it and saves as .png
save(savefilename,'Cp','t','AIFparams','FA','TR','bloodsignal')