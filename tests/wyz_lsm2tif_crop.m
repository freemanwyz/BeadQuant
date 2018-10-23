%Crop LSM and convert to multistack tiff
warning('off', 'Images:initSize:adjustingMag');
xrg = 301:500;
yrg = 301:500;

fnameBg = 'bgcrop.tif';
fnameTs = 'tscrop.tif';

%% load files
[fname,stub] = uigetfile('*.lsm','Select the background .lsm file');
path0 = fullfile(stub,fname);
[fname,stub] = uigetfile(fullfile(stub,'*.lsm'),'Select the .lsm file you are analyzing');
path1 = fullfile(stub,fname);

%% background
dat0 = tiffread(path0);
maxVal = 2^dat0.bits-1;
bgch1 = dat0.data{1};
imwrite(bgch1(xrg,yrg), fnameBg);
bgch2 = dat0.data{2};
imwrite(bgch2(xrg,yrg), fnameBg, 'WriteMode', 'append');

%% time lapse
dat1 = tiffread(path1);
nFrame = length(dat1);
sigch2 = zeros(Nx,Ny,nFrame);
tt = 1;
for ii=1:nFrame
    tmp = dat1(ii);
    datEle = tmp.data{1};
    if tt==1
        imwrite(datEle(xrg,yrg), fnameTs);
    else
        imwrite(datEle(xrg,yrg), fnameTs, 'WriteMode', 'append');
    end
    tt = tt + 1;
    datEle = tmp.data{2};
    imwrite(datEle(xrg,yrg), fnameTs, 'WriteMode', 'append');
    tt = tt + 1;
end