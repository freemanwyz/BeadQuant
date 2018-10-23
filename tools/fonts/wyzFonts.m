%% convert font images to digits
% dat = imread('inconsolate_digits.png');
dat = imread('calibri_digits.png');
dat = 255 - rgb2gray(dat);
dat(dat<50) = 0;
[Nx,Ny] = size(dat);
dat = [dat,zeros(Nx,100)];
wid = 91;
res = zeros(Nx,wid,10);
for ii=1:10
    ix = 11 + (ii-1)*wid;
    res(:,:,ii) = dat(:,ix:(ix+wid-1));
    figure;imshow(res(:,:,ii))
    pause(0.2)
    close
end

%% resize
rzRatio = 0.25;
tmp = imresize(res(:,:,1),rzRatio);
[Nx5,Ny5] = size(tmp);
res5a = zeros(Nx5,Ny5,10);
xrg = zeros(10,2);
for ii=1:10
    ix = 11 + (ii-1)*wid;
    res5a(:,:,ii) = imresize(dat(:,ix:(ix+wid-1)),rzRatio);
    
    [px,py] = find(res5a(:,:,ii)>0);
    xrg(ii,1) = min(px);
    xrg(ii,2) = max(px);    
    
    figure;imshow(res5a(:,:,ii))
    pause(0.2)
    close
end

xmin = min(xrg(:,1));
xmax = max(xrg(:,2));
res5 = zeros(xmax-xmin+1,Ny5,10);
for ii=1:10
    res5(:,:,ii) = res5a(xmin:xmax,:,ii);
    figure;imshow(res5(:,:,ii))
    pause(0.2)
    close
end

save('./tools/fonts/digitTemp.mat','res','res5');
