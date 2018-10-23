%% Bead detection based on circular Hough Transform
bgMean = dd.bgMean;  % Please first load the session file from beadGUI for this option
cropImg = 0;
addNoise = 0;
if cropImg
    imga = bgMean(300:500,300:500);
else
    imga = bgMean;
end
if addNoise
    image_intensity0 = 0:0.01:1;
    varMax = 0.3;
    var0 = 0:varMax/100:varMax;
    imgb = imnoise(imga,'localvar',image_intensity0,var0);
else
    imgb = imga;
end
d0a = sqrt(imgb);
% d0a = imgb;
% d0a = imgaussfilt(d0a,3);
% d0a = wiener2(d0a,[7,7]);
% d0a = wiener2(d0a,[5,5]);

%% bead detection using Hough transform
xrg = [8 20];
% xrg = [3,15];
[centers, radii, metric, accumMatrix] = imfindcircles0(d0a,xrg);
% [centers, radii, metric] = imfindcircles(d0a,xrg,'ObjectPolarity','bright','Sensitivity',0.95, 'EdgeThreshold',0.02);
% K0 = cat(3,imgb*0,imgb,imgb*0);
% figure;imshow(K0);
% viscircles(centers, radii,'EdgeColor','b');

%% plot
% K0 = cat(3,bgMean*0,d0a,imgb*0);
K0 = cat(3,bgMean*0,imgb,imgb*0);
col = {'blue','red',[0.99 0.99 0.99],'yellow',[1 0.5 0],[0.99 0 0.99]};
for ii=1:length(centers)
    centersAll = centers{ii};
    radiiAll = radii{ii};
    K0 = insertShape(K0,'circle',[centersAll,radiiAll],'Color',col{ii},'Opacity',1,'SmoothEdges',false);
end
imshow(K0);

K1 = cat(3,bgMean*0,imgb,imgb*0);
for ii=1:length(centers)
    centersAll = centers{ii};
    radiiAll = radii{ii};
    K1 = insertShape(K1,'circle',[centersAll,radiiAll],'Color',[0.9 0.9 0.9],'Opacity',0.4,'SmoothEdges',false);
end
imwrite(K1,'bead_073115_fast_0.05_wiener_0.tiff');
% imwrite(K1,'bead_060815_fast.tiff');

% a0 = abs(accumMatrix{4});
% K0 = cat(3,a0/max(a0(:)),imga*0.5,imga*0);
% figure;imshow(K0);



