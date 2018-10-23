function varargout = beadGUI(varargin)
%BEADGUI MATLAB code for beadGUI.fig

% Last Modified by GUIDE v2.5 22-Jul-2015 10:56:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @beadGUI_OpeningFcn, ...
    'gui_OutputFcn',  @beadGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% init -----
function beadGUI_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
addpath(genpath('./src/'))
handles = beadAddPath(handles);
movegui(gcf,'center')
guidata(hObject, handles);
% plot something in the axes 
% tmp = rand(100,100);

% [X,map] = imread('./cfg/beads.png');
% if ~isempty(map)
%     Im = ind2rgb(X,map);
% end
% imshow(Im);

Im = imread('./cfg/beads.png');
imshow(Im);
[Nx,Ny,~] = size(Im);
text(round(Nx/5),round(Ny/2.5),'BEAD','FontSize',128,'Color',[1 0.5 0])
text(round(Nx/5)+15,round(Ny/1.7),'Version 0.72','FontSize',32,'Color',[1 0.5 0])

% detection -----
function handles = detectBeads_Callback(hObject, ~, handles)
if isempty(handles.dd.img0) && isempty(handles.dd.img1)
    msgbox('Please select image first');
else
    set(handles.done1,'String', 'Busy');
    pause(0.1);
    handles = beadGatherParam(handles);
    t0 = tic;
    handles = beadDetectAlignExtract( handles );
    t1 = toc(t0);
    fprintf('This step takes %f seconds\n',t1);
    handles.pp.doneStage = 1;
    set(handles.done1,'String', 'Done');    
    guidata(hObject, handles);
end

function handles = fitCurve_Callback(hObject, ~, handles)
if handles.pp.doneStage >=1
    set(handles.done2,'String', 'Busy');
    pause(0.1);
    handles = beadGatherParam(handles);
%     handles = beadProcessor0(handles);
    handles = beadProcessor(handles);
    set(handles.done2,'String', 'Done');
    handles.pp.doneStage = 2;
    guidata(hObject, handles);
else
    msgbox('Please detect beads first','Info');
end

function handles = chooseBeads_Callback(hObject, ~, handles)
if handles.pp.doneStage >=2
    set(handles.done3,'String', 'Busy');
    pause(0.1);
    handles = beadGatherParam(handles);
    dd = handles.dd;
    pp = handles.pp;
    tic
    controlBeadIdx = beadControlIdx(dd.resBead,dd.bgMean,dd.controlBeadCenters);
    goodBeadIdx = beadChoose(dd.dff_fit,dd.tau_fit,dd.dff_mat,dd.time_vect,...
        controlBeadIdx,pp.useDff,pp.fitModel,pp.outPath,handles);
    beadsDraw(dd.resBead,dd.bgMean,dd.sig,goodBeadIdx,...
        controlBeadIdx,pp.outPath,pp.maxVal)
    handles.dd.goodBeadIdx = goodBeadIdx;
    handles.dd.controlBeadIdx = controlBeadIdx;
    beadInfoWrite(handles.dd, pp.outPath);
    toc
    set(handles.done3,'String', 'Done');
    handles.pp.doneStage = 3;
    guidata(hObject, handles);
    set(handles.myinfo,'String','All Done!');
else
    msgbox('Please detect beads and fit curve first','Info');
end

function runAll_Callback(hObject, eventdata, handles)
handles = detectBeads_Callback(hObject, eventdata, handles);
set(handles.done1,'String', 'Done');
handles = fitCurve_Callback(hObject, eventdata, handles);
set(handles.done2,'String', 'Done');
handles = chooseBeads_Callback(hObject, eventdata, handles);
set(handles.done3,'String', 'Done');
handles.pp.doneStage = 3;
guidata(hObject, handles);

% load files and set output path -----
function open1_Callback(hObject, ~, handles)
handles = beadOpenImg(handles,0);
set(handles.myinfo,'String','One image loaded, load another or continue');
guidata(hObject, handles);

function open2_Callback(hObject, ~, handles)
handles = beadOpenImg(handles,1);
set(handles.myinfo,'String','Two images loaded, please continue');
guidata(hObject, handles);

function outFolder_Callback(hObject, ~, handles)
folder_name = uigetdir('output','Ouput folder');
handles.pp.outPath = folder_name;
if length(folder_name)>40
    folder_name = ['..',folder_name(end-20:end)];
end
set(handles.outFolderName,'String',folder_name);
guidata(hObject, handles);

% tools -----
function loadSession_Callback(hObject, ~, handles)
[fname,stub] = uigetfile('*.session.mat','Select session file');
if length(stub)>1
    h0 = msgbox('Load session...');
    pathSession = fullfile(stub,fname);
    tmp = load(pathSession);
%     tmp.dd.img0 = double(tmp.dd.img0)/65535;
%     tmp.dd.img1 = double(tmp.dd.img1)/65535;
%     tmp.dd.sig = double(tmp.dd.sig)/65535;
    if ~isfield(tmp.pp,'smoMethod')  % compatibility
        tmp.pp.smoMethod = 1;
        tmp.pp.thrDetect = 0.05;
        tmp.pp.thrMulti = 3;
        tmp.pp.maxVal = 2^16-1;
    end
    if ~isfield(tmp.pp,'extractMethod')
        tmp.pp.extractMethod = 1;
    end
    handles.pp = tmp.pp;
    handles.dd = tmp.dd;
    beadSetParam(handles)
    guidata(hObject, handles);
    close(h0);
end

function saveSession_Callback(~, ~, handles)
[file,path] = uiputfile('*.*','Choose folder to save the session and enter file name');
if length(path)>1
    h0 = msgbox('Saving session...');
    handles = beadGatherParam(handles);
    pp = handles.pp;
    dd = handles.dd;
    if isa(dd.sig,'double')  % compatibility
        dd.img0 = uint16(dd.img0*65535);
        dd.img1 = uint16(dd.img1*65535);
        dd.sig = uint16(dd.sig*65535);
    end
    save([path filesep file '.session.mat'],'pp','dd','-v7.3');
    close(h0);
end

function startParpool_Callback(~, ~, ~)
ver0 = version('-release');
ver0 = str2double(ver0(1:4));
if ver0>=2013
    p0 = gcp('nocreate');
    if isempty(p0)
        parpool;
    else
        msgbox('Parallel pool already running');
    end
else
    msgbox('Please start parallel pool manually');
end

function varargout = beadGUI_OutputFcn(~, ~, handles)
varargout{1} = handles.output;


