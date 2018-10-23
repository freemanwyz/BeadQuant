%% Joey
if 0
    fname = 'D:\neuro_WORK\coculture_ucdavis_WORK\MFD_dat\040115\0.02ColNoPEG\20xTS_100uMGluASide@T=50.lsm';
    % fname = 'D:\neuro_WORK\coculture_ucdavis_WORK\MFD_dat\040715\Col0.05NoPEG\20xTS_SponN&ASameSide.lsm';
    lsminfo(fname)
end

%% GECI
ftop = 'D:\neuro_WORK\coculture_ucdavis_WORK\';

ii = 5;
fprintf('%d ----------\n',ii);
fname = [ftop 'GECI\' 'nn.4001-sv' num2str(ii) '.avi'];
fnameOut = [ftop 'GECI_tiff\' 'nn.4001-sv' num2str(ii) '.tif'];
xyloObj = VideoReader(fname);
vidWidth = xyloObj.Width;
vidHeight = xyloObj.Height;
vidLen = xyloObj.Duration*xyloObj.FrameRate;

% while hasFrame(xyloObj)
for kk = 1:vidLen
    fprintf('%d\n',kk);
    x0 = rgb2gray(readFrame(xyloObj));
    imwrite(x0,fnameOut,'WriteMode', 'append',  'Compression', 'none');
end

if 0
    mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
        'colormap',[]);
    k = 1;
    while hasFrame(xyloObj)
        mov(k).cdata = readFrame(xyloObj);
        k = k+1;
    end
    hf = figure;
    set(hf,'position',[150 150 vidWidth vidHeight]);
    movie(hf,mov,1,xyloObj.FrameRate);
end


