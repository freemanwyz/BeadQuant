function [ flst, flstName ] = getFlst( file )
%getFlst get file list
%   fpath: path/name, fname: name only
fid = fopen(file);
flst = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
flst = flst{1};
flstName = cell(20,1);
for ii=1:length(flst)
    t1 = strsplit(flst{ii},'/');
    flstName{ii} = t1(end);
end
flst = cellfun(@(x) strrep(x,'/',filesep), flst,'UniformOutput',0);
% flstName = cellfun(@(x) strrep(x,'/',filesep), flstName,'UniformOutput',0);
end

