function [bgch1_rigid,bgch2_rigid,qual] = beadAlignment(sigch1Mean, bgch1, bgch2, outPath, mthd)
% BEAD_REGISTRATION_CURVE Align images and extract curves

[Nx,Ny] = size(bgch1);
rt0 = mean(bgch1(:))/mean(sigch1Mean(:));
bgch1 = bgch1/rt0;

% % some preprocessing
% xMid = round(Nx/2);
% yMid = round(Ny/2);
% sz = 1000;
% 
% rgx = (xMid-sz):(xMid+sz);
% rgy = (yMid-sz):(yMid+sz);
% %     rgx = 1:2*sz;
% %     rgy = 1:2*sz;
% rgx = rgx(rgx>0 & rgx<=Nx);
% rgy = rgy(rgy>0 & rgy<=Ny);
% bgch1Crop = bgch1(rgx ,rgy );
% 
% % bgch1Crop = wiener2(bgch1Crop,[5 5]);
% sigch1MeanCrop = sigch1Mean( rgx ,rgy );

USE_RESIZE = 1;
if USE_RESIZE
    if Nx>6000
        RESIZE_RATIO = 0.125;
    elseif Nx>2500
        RESIZE_RATIO = 0.25;
    else
        RESIZE_RATIO = 0.5;
    end
    bgch1Crop = imresize(bgch1,RESIZE_RATIO);
    sigch1MeanCrop = imresize(sigch1Mean,RESIZE_RATIO);
    bgch1Crop(bgch1Crop<0) = 0;
    sigch1MeanCrop(sigch1MeanCrop<0) = 0;
    bgch1Crop = wiener2(bgch1Crop,[5 5]);
    sigch1MeanCrop = wiener2(sigch1MeanCrop,[5 5]);
    %     bgch1Crop = sqrt(bgch1Crop);
    %     sigch1MeanCrop = sqrt(sigch1MeanCrop);
end

%% IAT, slow
if mthd==1
    tic;
    par.transform = 'euclidean';
    par.levels = 2;
%     par.iterations = 20; %iterations per level
    d1Warp = iat_LucasKanade(bgch1Crop, sigch1MeanCrop, par);
%     d1Warp = iat_LucasKanade(bgch1, sigch1Mean, par);
    if USE_RESIZE
        d1Warp(:,3) = d1Warp(:,3)/RESIZE_RATIO;
    end

    bgch1_rigid = iat_inverse_warping(bgch1, d1Warp, par.transform, 1:Ny, 1:Nx);
    Ka = cat(3,bgch1_rigid,sigch1Mean,bgch1_rigid*0);
    imwrite(double(Ka),[outPath,filesep,'alignment_iat.tiff']);
    [bgch2_rigid, ~] = iat_inverse_warping(bgch2, d1Warp, par.transform, 1:Ny, 1:Nx);
    t1 = toc;
    fprintf('Time for alignment is %f second\n',t1);
end

%% Matlab
if mthd==0    
    tic;
    [optimizer, metric]  = imregconfig('multimodal');
%     optimizer.MaximumIterations = 500;
%     optimizer.RelaxationFactor = 0.75;
    optimizer.GrowthFactor = 1.01;
%     optimizer.InitialRadius = 0.005;
    optimizer.Epsilon = 1e-8;
    optimizer.MaximumIterations = 1000;
    
%     tr0 = imregtform(bgch1,sigch1Mean,'rigid',optimizer, metric);
    tr0 = imregtform(bgch1Crop,sigch1MeanCrop,'rigid',optimizer, metric);   
    tr0.T
    if USE_RESIZE
        t0 = tr0.T;
        t0(3,1:2) = t0(3,1:2)/RESIZE_RATIO;
        tr0.T = t0;
    end
    
    bgch1_rigid = imwarp(bgch1,tr0,'OutputView',imref2d(size(sigch1Mean)));
    bgch2_rigid = imwarp(bgch2,tr0,'OutputView',imref2d(size(sigch1Mean)));
    
    Ka = cat(3,bgch1_rigid,sigch1Mean,bgch1_rigid*0);
    imwrite(double(Ka),[outPath,filesep,'alignment_matlab.tiff']);
    
    t1 = toc;
    fprintf('%f\n',t1);
    
    % d1Registeredx = imwarp(d1,d1tform);
    % d1Registered = imregister(d1,d0,'rigid',optimizer, metric);
    
    % figure;imshowpair(d0, d1,'Scaling','joint');
    % figure;imshowpair(d0, d1_trans,'Scaling','joint');
    % figure;imshowpair(d0, d1_rigid,'Scaling','joint');    
end

x1 = bgch1_rigid(:);
x2 = sigch1Mean(:);
% c0 = corrcoef(x1,x2);
% qual = c0(1,2);
qual = mutualinfo(double(x1),double(x2));

fprintf('Quality %f\n',qual);


