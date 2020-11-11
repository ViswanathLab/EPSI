%% Script for opening Bruker EPSI data
% Input: EPSI acquisition scan number (data_path) and subfolder (path_img)
% for Magnitude data (1) or Phase (2)
% Output: MetImage - EPSI data as 4D matrix
% GB - 20190401

data_path =  './8';     % put the data folder here
path_img   = [data_path '/pdata/1/'];
path_method  = [data_path '/method'];
path_recon = [path_img 'reco']; 
%% read recon header
recon_header = readBrukerHeader(path_recon);   
%% read method header
method_header = readBrukerHeader(path_method);
n_rep = method_header.PVM_NRepetitions;%6;  
dimension = [recon_header.RECO_size; n_rep]'; % [x spec y N_d] the dimension of the 2dseq file
%% read 2dseq file
f1 = fopen([path_img '2dseq'],'r');
im_all=fread(f1,'int16');
im_epsi=reshape(im_all,dimension(1),dimension(2),dimension(3),dimension(4));
%% rotate epsi image
resize_fac = 1;%resize factor
im_rot = permute(im_epsi(:,:,:,:),[3 1 2 4]);
im_rot = flip(im_rot,1);
inim = imresize(sum(sum(im_rot,3),4),[dimension(1) dimension(3)]*resize_fac);
% figure, imshow(inim,[]), title('image summed over all dynamics and spectral points')
%% compute the ppm range for spectra
cent_ppm = method_header.PVM_FrqWorkOffsetPpm;
rang_ppm = method_header.SpecBandPpm;
npts     = dimension(2);
start_ppm = cent_ppm + rang_ppm/2;
stop_ppm  = cent_ppm - rang_ppm/2;
step_ppm  = rang_ppm/(npts-1);

idxarry_ppm = -start_ppm:step_ppm:-stop_ppm;
idxarry_ppm = -idxarry_ppm;
%% RUN until here
MetImage = im_rot;      % if you have folder 1 (magnitude)
%%
MetImageP = im_rot;     % if you have folder >1 (phase)
MetImageC = MetImage.*exp(1i*MetImageP);
