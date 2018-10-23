function p = getInitParamBead( varargin )
%GETINITPARAM Set the parameter to synpase detection using profiles
%   preset: typical, simFP, simPR

ipa = inputParser;
addOptional(ipa,'preset','typical');
parse(ipa,varargin{:})
up = ipa.Results;

p = [];

% get path -----
if ispc
  pathHome = [getenv('HOMEDRIVE') getenv('HOMEPATH')]; 
else 
  pathHome = getenv('HOME'); 
end
myCfgFile = [pathHome filesep 'wyz_dat'];
fileID = fopen(myCfgFile,'r');
d0 = textscan(fileID,'%s','Delimiter','=');
idx = find(cellfun(@(x) strcmp(x,'ucdavis_grace_bead'), d0{1}));
if isempty(idx)
    warning('Can not find data and dump file path\n');
    keyboard
end
dpath = d0{1}{idx+1};
dpath = strrep(dpath,'/',filesep);
fclose(fileID);
p.tp = dpath;  % top path
p.SE = strel('square',3);

% presets -----
if strcmp(up.preset,'typical')
    p.crop = 1;
    p.crop_rgx = 1801:2400;
    p.crop_rgy = 1701:2200;
end


end

