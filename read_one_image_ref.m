%% Read and save the reference anatomical HR image
% Input: HR acquisition scan number (data_path)
% Output: RefImage - HR image
% GB - 20190401

clear all
close all

data_path =  './2';     % HR anatomical sequence
path_img   = [data_path '/pdata/1/']; 
path_method  = [data_path '/method'];
path_recon = [path_img 'reco'];               
%%
recon_header = readBrukerHeader(path_recon);   
dimension = [recon_header.RECO_size; recon_header.RecoObjectsPerRepetition]'; % [x y z] the dimension of the 2dseq file
%% read method header
method_header = readBrukerHeader(path_method);
n_rep = method_header.PVM_NRepetitions;%6;                
%% read 2dseq file
f1 = fopen([path_img '2dseq'],'r');
im_all=fread(f1,'int16');
%%
image = reshape(im_all, [dimension, n_rep]);
%%
for i = 1:size(image,3)
    im_rot = flip(permute(sum(squeeze(image(:,:,i)),3),[2 1]),1);
    %im_rot = imgaussfilt(im_rot,1);
    figure(i), imshow(im_rot,[])
end
%% Save salected image
RefImage = flip(permute(sum(squeeze(image(:,:,9)),3),[2 1]),1);
figure, imshow(RefImage,[])