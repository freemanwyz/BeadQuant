function handles = beadAddPath(handles)

addpath(genpath('./cfg/'))
% addpath(genpath('./examples/'))
addpath(genpath('./output/'))
addpath(genpath('./tools/'))
% Java files for writing XLS file
javaaddpath('tools/xlwrite/poi_library/poi-3.8-20120326.jar');
javaaddpath('tools/xlwrite/poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('tools/xlwrite/poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('tools/xlwrite/poi_library/xmlbeans-2.3.0.jar');
javaaddpath('tools/xlwrite/poi_library/dom4j-1.6.1.jar');
javaaddpath('tools/xlwrite/poi_library/stax-api-1.0.1.jar');
warning('off', 'Images:initSize:adjustingMag');
if exist('output','dir')~=7
    mkdir('output');
end

% read default parameters
tmp = load('./cfg/cfgDefault.mat');
handles.pp = tmp.pp;
handles.dd = tmp.dd;

end