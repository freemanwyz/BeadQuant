% TODO: support more input image numbers
function handles = beadOpenImg(handles,fcnt)
% fspec = {'*.lsm';'*.tif';'*.tiff'};
[fname,stub] = uigetfile('*.*','Select File');

if isempty(stub) || length(stub)<2
    return
end

tmp = regexp(fname,'[.]','split');
imgType = tmp{end};
fullName = fullfile(stub,fname);

% set image parameter for LSM
if strcmp(imgType,'lsm')
    [lsmA,~,~] = lsminfo(fullName);
    nFrame = lsmA.TIMESTACKSIZE;
    nChan = lsmA.NUMBER_OF_CHANNELS;
    r0 = round(44e-6/lsmA.VoxelSizeX);
    r1 = round(154e-6/lsmA.VoxelSizeX);
    if r0<3
        r0 = 3;
    end
    if r1<r0
        r1 = r0;
    end
    set(handles.radRg0,'String',num2str(r0));
    set(handles.radRg1,'String',num2str(r1));
    if nFrame>1
        handles.pp.timeStep = round(lsmA.TimeStamps.AvgStep);
        set(handles.timeGap,'String',num2str(handles.pp.timeStep));
    end
    h = waitbar(0,'Reading Image...');
    dat1 = tiffread(fullName);
    waitbar( 10 / (nFrame+10))
    maxVal = 2^dat1(1).bits-1;
    tmp = dat1(1);
    datShow = tmp.data{1};
    [Nx,Ny] = size(tmp.data{1});
    imgx = zeros(Nx,Ny,nFrame,nChan,'uint16');
%     imgx = zeros(Nx,Ny,nFrame,nChan);
    for ii=1:nFrame
        tmp = dat1(ii);
        for jj=1:nChan
%             imgx(:,:,ii,jj) = double(tmp.data{jj})/maxVal;
%             imgx(:,:,ii,jj) = double(tmp.data{jj})/maxVal;
            imgx(:,:,ii,jj) = tmp.data{jj};
            imgx(:,:,ii,jj) = tmp.data{jj};
        end
        h = waitbar( (ii+10) / (nFrame+10));
    end
    if isvalid(h)
        close(h)
    end
else  % read TIFF after user select channel number
    tmp = imfinfo(fullName);
    maxVal = 2^tmp(1).BitDepth-1;
    nAll = length(tmp);
    if isempty(handles.pp.nChan)
        nChan = str2double(inputdlg('Please set channel number:','Input'));
    else
        nChan = handles.pp.nChan;
    end
    nFrame = nAll/nChan;
    [Nx,Ny] = size(imread(fullName,1));
    datShow = imread(fullName,1);
    imgx = zeros(Nx,Ny,nFrame,nChan,'uint16');
%     imgx = zeros(Nx,Ny,nFrame,nChan);
    h = msgbox('Reading image');
    for ii=1:nFrame
        for jj=1:nChan
            imgx(:,:,ii,jj) = imread(fullName,(ii-1)*nChan+jj);
%             imgx(:,:,ii,jj) = double(imread(fullName,(ii-1)*nChan+jj))/maxVal;
        end
    end
    if isvalid(h)
        close(h)
    end
    if fcnt==1
        msgbox('Please set the time step and radius range');
        uiwait
    end
end

imshow(datShow);

% update parammeters and data
handles.pp.nChan = nChan;
handles.pp.maxVal = maxVal;
set(handles.channelNum,'String',num2str(nChan));

str0 = {{'1'},{'1','2'},{'1','2','3'}};
set(handles.selbgch,'String',str0{nChan});
set(handles.selsigch,'String',str0{nChan});

if nFrame>1    
    handles.pp.timeRg0 = 1;
    handles.pp.timeRg1 = nFrame;
    set(handles.time1,'String',num2str(handles.pp.timeRg1));
end

if fcnt==0
    handles.dd.img0 = imgx;
else
    handles.dd.img1 = imgx;
end

% handles.pp.stub = stub;
if fcnt==0
    handles.pp.ftype0 = imgType;
    handles.pp.path0 = fullName;
    handles.pp.fname0 = fname;
    set(handles.img0Name,'String',handles.pp.fname0);
else
    handles.pp.ftype1 = imgType;
    handles.pp.path1 = fullName;
    handles.pp.fname1 = fname;
    set(handles.img1Name,'String',handles.pp.fname1);
end

end


