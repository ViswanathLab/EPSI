%% code to process EPSI data
% - zero fill & apodize spectra
% - plot spectra: every time frame all voxels, single voxel, overlay
% different voxels
% - create metabolic maps by integrating peaks and overlay them on
% anatomical image
% evaluate SNR of defined ROI and ratio of metabolites
% GB - 20190917

%% zero fill + apodization
% fill the number of points for zero filling and apodization width
ZeroFillPoints = 0;     % values: 0 or 512 or 1024
Apod = 0;               % values 0 to 15

yiC_FID = fft(MetImageC(:,:,:,:),size(MetImageC,3),3);
yiC_FID_z = zeros(8,8,ZeroFillPoints+size(MetImageC,3),size(MetImageC,4));
yiC_FID_z(1:8,1:8,1:size(MetImageC,3),1:size(MetImageC,4)) = yiC_FID;

time_z = [0:(size(yiC_FID_z,3)-1)]/(19.47*32);

smoothFilt = exp(-time_z*Apod*pi)*10;      % apodization
smoothFilt4D = repmat(smoothFilt',[1 8 8 size(MetImageC,4)]);
smoothFilt4D = permute(smoothFilt4D,[2,3,1,4]);

yiC_FID_zf = yiC_FID_z.*smoothFilt4D;
yiC_zf = fft(yiC_FID_zf,ZeroFillPoints+size(MetImageC,3),3);

idxarry_ppm_zf = linspace(min(idxarry_ppm),max(idxarry_ppm),size(yiC_zf,3));

%% place points on HR image for voxel reference (grid on anatomical)
figure
imshow(RefImage,[])
hold on
plot(ones(8,1)*[1:8]*32,[1:8]'*ones(1,8)*32,'.r')

%% Spectra visualization
%% Plot spectra over 10 time frames
xy = [1;8;1;8];
for tr = 1:10
figure(tr)
n = 1;
mApl = max(max(max(sum(abs(yiC_zf(xy(1):xy(2),xy(3):xy(4),:,tr)),4))));
for i = xy(1):xy(2)
    for j = xy(3):xy(4)
        subplot((-xy(1)+xy(2)+1),(xy(4)-xy(3)+1),n)
        plot(idxarry_ppm_zf,squeeze(sum(abs(yiC_zf(i,j,:,tr)),4)))
        set(gca, 'XDir','reverse');
        ylim([0 mApl])
%         xlim([170 188])
        n=n+1;
        hold off
    end
end
end

% %% Plot one spectra for one time point (line, column, ppm, time frame)
% figure, plot(idxarry_ppm_zf,squeeze(abs(yiC_zf(3,3,:,4)))) % if u want to plot more time frames use (start : stop)
% 
% %% if you want to plot controlateral brain spectra and tumor
% figure, plot(idxarry_ppm_zf,squeeze(abs(yiC_zf(3,3,:,1:10))./max(abs(MetImageC(3,3,:,1:10)))))
% hold on
% plot(idxarry_ppm_zf,squeeze(abs(yiC_zf(3,4,:,1:10))./maxabs(MetImageC(3,4,:,1:10))))
% %% plot one spectra as sum of time points
% figure, plot(idxarry_ppm_zf,squeeze(sum(abs(yiC_zf(3,4,:,2:10)),4)))   % line, column, ppm, time frames (start : end)
% 
% %% plot one spectra as sum of voxels in column, different lines
% figure, plot(idxarry_ppm_zf,squeeze(sum(abs(yiC_zf(3:4,3,:,2)),1)))   % lines (start:end), column, ppm, time frame
% 
% %% plot one spectra as sum of voxels in line, different column
% figure, plot(idxarry_ppm_zf,squeeze(sum(abs(yiC_zf(3,[3:4],:,2)),2)))   % line, columns (start:end), ppm, time frame
% 
% %% plot one spectra as sum of any voxels with any voxel - line, columns, ppm, time frame
% figure, plot(idxarry_ppm_zf,squeeze(abs(yiC_zf(3,4,:,2) + yiC_zf(3,4,:,3)))) 

%% Metabolic maps
%% Pyruvate - integration of peak
for i = 1:10
ppm_a = 168; %ppm start point for integration 
ppm_b = 173; %ppm stop point for integration
spectra_range = find((idxarry_ppm_zf>ppm_a)&(idxarry_ppm_zf < ppm_b));
dyn_range = i;%1:n_rep;% dynamic range for integration

inim_pyr(:,:,i) = imresize(sum(sum(abs(yiC_zf(:,:,spectra_range,dyn_range)),3),4),[dimension(1) dimension(3)]*resize_fac);
end

%% Lactate - integration of peak
for i = 1:10
ppm_a = 182; %ppm start point for integration 
ppm_b = 184; %ppm stop point for integration
spectra_range = find((idxarry_ppm_zf>ppm_a)&(idxarry_ppm_zf < ppm_b));
dyn_range = i;%1:n_rep;% dynamic range for integration

inim_lac(:,:,i) = imresize(sum(sum(abs(yiC_zf(:,:,spectra_range,dyn_range)),3),4),[dimension(1) dimension(3)]*resize_fac);
end

%% Alanine - integration of peak
for i = 1:10
ppm_a = 174; %ppm start point for integration 
ppm_b = 178; %ppm stop point for integration
spectra_range = find((idxarry_ppm_zf>ppm_a)&(idxarry_ppm_zf < ppm_b));
dyn_range = i;%1:n_rep;% dynamic range for integration

inim_ala(:,:,i) = imresize(sum(sum(abs(yiC_zf(:,:,spectra_range,dyn_range)),3),4),[dimension(1) dimension(3)]*resize_fac);
end

%% create mask for brain-region of interest
figure
imshow(RefImage,[])
set(gcf,'position',[10,10,1000,1000])
message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
uiwait(msgbox(message));
hFH = imfreehand(); % Actual line of code to do the drawing.
% Create a binary image ("mask") from the ROI object.
RefMask = hFH.createMask();

%% overlay one metabolic image as ratio of metabolites on anatomical
figure
B = repmat(mat2gray(double(RefImage)),[1,1,3]);

% Display the back image
hB = imagesc(B);axis image off;

% Add the front image on top of the back image
hold on;
F = imresize(inim_lac(:,:,7),[256 256],'lanczos2')./imresize(inim_pyr(:,:,7),[256 256],'lanczos2');
  % change time frame 5 to the one needed

F_Mask = F.*RefMask;
hF = imagesc(F_Mask,[min(F_Mask(:)) max(F_Mask(:))]);

% Make the foreground image transparent
alphadata = 0.5.*(F >= min(F(:)));
set(hF,'AlphaData',alphadata);
% colorbar

%% Maps of metabolite to STD of noise
%% Evaluate noise (as the STD of the real part of spectra of voxel (1,1))
for i = 1:20
    inim_noise(:,:,i) = std(real(yiC_zf(:,:,:,i)),0,3);
end

%% overlay metabolic image on anatomical for 9 frames (Ratio to STD value of noise)
figure
I = repmat(mat2gray(double(RefImage)),[1,1,3]);

% metabolite to be plot (change inim_? to metabolite, ?=lac or glc)
MetPlot = inim_pyr;
Fnoise = squeeze(inim_noise(1,1,:));

for i=1:9       % number of time frames - do not forget to change dimantion in subplot below

subplot(3,3,i)

F = imresize(MetPlot(:,:,i),[256 256],'lanczos2')./Fnoise(i);
F_Mask = F.*RefMask;

% Display the back image
hB = imagesc(I);axis image off;

% Add the front image on top of the back image
hold on;
F_MinMax = [min(min(F_Mask)) max(max(F_Mask))];
hF = imagesc(F_Mask,F_MinMax);

% Make the foreground image transparent
alphadata = 0.6.*(F_Mask >= F_MinMax(1));
set(hF,'AlphaData',alphadata);
colorbar
end

%% overlay one metabolic image as ratio of metabolite to std of noise on anatomical
figure
B = repmat(mat2gray(double(RefImage)),[1,1,3]);

% Display the back image
hB = imagesc(B);axis image off;

% Add the front image on top of the back image
hold on;

% metabolite to be plot (change inim_? to metabolite, ?=lac or glc)
MetPlot = inim_lac(:,:,5);
Fnoise = inim_noise(1,1,5);
F = imresize(MetPlot,[256 256],'lanczos2')./Fnoise;


F_Mask = F.*RefMask;
hF = imagesc(F_Mask,[min(F_Mask(:)) max(F_Mask(:))]);

% Make the foreground image transparent
alphadata = 0.5.*(F >= min(F(:)));
set(hF,'AlphaData',alphadata);
% colorbar

%% Define voxel to calculate metabolite levels
figure, imagesc(RefImage)
set(gcf,'position',[10,10,500,500])
message = sprintf('Move the voxel at the area of interest (tumor/control)');
uiwait(msgbox(message));
hFH = imrect(gca,[10 10 15 15]);
position = wait(hFH);
RefMaskMet = hFH.createMask();

%% calculate ratio of product to substrate (select frame (frame...) and compounds (inim_...))
frame_Sub = 4;
frame_Prod = 4;

inim_Sub = inim_ala;
inim_Prod = inim_lac;

vProd = mean(mean(imresize(inim_Prod(:,:,frame_Prod),[256 256],'lanczos2').*RefMaskMet));
vSub = mean(mean(imresize(inim_Sub(:,:,frame_Sub),[256 256],'lanczos2').*RefMaskMet));
ratio = vProd./vSub

%% calculate ratio of metabolite to noise (select frame (frame...) and compounds (inim_...))
frame_Met = 4;

inim_Met = inim_lac;

inim_SNR = imresize(inim_Met(:,:,frame_Met),[256 256],'lanczos2')/inim_noise(1,1,frame_Met);

SNR = sum(sum(inim_SNR.*RefMaskMet))/15/15

