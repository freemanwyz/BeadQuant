%% put soma to synD
% use 8 bit three channel LSM as input
% freemanwyz@gmail.com
warning('off', 'Images:initSize:adjustingMag');
clear all;
close all;

thr0 = 40;
min_size = 100;

fid = fopen('f1.txt');
flst = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
flst = flst{1};

ii = 17;

%% detection
for ii=1:length(flst)
    % read image
    fin = ['./',flst{ii},'.lsm'];
    fprintf('%s  >>>>>>>>>>\n',fin);
    [I, Iarg] = lsmread(fin);
    Nx = Iarg.dimX;
    Ny = Iarg.dimY;
    Ineu = I(1,1,1,:,:);
    Isyn = I(1,2,1,:,:);
    Ituj = I(1,3,1,:,:);
    Ineu = double(reshape(Ineu,[Nx,Ny]))/256;
    Isyn = double(reshape(Isyn,[Nx,Ny]))/256;
    Ituj = double(reshape(Ituj,[Nx,Ny]))/256;
    Ituj(Ituj<thr0/256) = 0;
    I1 = cat(3,Ineu,Isyn,Ituj);
    fname_out = ['./',flst{ii},'_tujthr.tif'];
    imwrite(I1,fname_out);
end

%%
% figure;imshow(Isyn);
% figure;imshow(Ituj);
% figure;imshow(Ituj>0.2);



