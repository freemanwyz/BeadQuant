function [ Img, ImgCrop, q ] = loadData( fpath, filetype, p )
%LOADDATA Load LSM image

q = [];
if strcmp(filetype,'lsm')
%     fin = [p.tp,filesep,'data',filesep,fpath,'.lsm'];
    fin = fpath;
    [I, Iarg] = lsmread(fin);
    q.Lx = Iarg.voxSizeX*1e6;
    Nx = Iarg.dimX;
    Ny = Iarg.dimY;
    Img = double(reshape(I(1,p.ch,1,:,:),[Nx,Ny]));    
else
    Img = imread(fpath);
    q.Lx = 0.2;
end

if p.crop
    ImgCrop = Img( p.crop_rgx, p.crop_rgy );
else
    ImgCrop = Img;
end
[Nx,Ny] = size(ImgCrop);
q.Nx = Nx;
q.Ny = Ny;
q.ntry = 5;
q.Lx2 = q.Lx^2;
q.Pix_per_synapse = round(1/q.Lx2);

Img = double(Img);
ImgCrop = double(ImgCrop);

end

