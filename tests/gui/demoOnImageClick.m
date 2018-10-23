clc;clear;
f0 = figure;
% hAxes = axes();
imObj = rand(500,500);
imageHandle = imshow(imObj);
sigma0 = 0.5;
set(imageHandle,'ButtonDownFcn',{@Click_CallBack sigma0});
set(f0,'CloseRequestFcn',@Close_CallBack);
% process_plot()
% set(imageHandle,'ButtonDownFcn',{@ImageClickCallback,sigma0,gca});

% FIGDisp = 1;%loop until user happy with rectangle
% while FIGDisp==1
%     imObj = rand(500,500);
%     f = figure;
%     imshow(imObj);
%     k = waitforbuttonpress;
%     point1 = get(gca,'CurrentPoint');    % button down detected
%     close(f);
% end