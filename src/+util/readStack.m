function [dat,datRef,dats] = readStack(path1,crop)

dat1 = tiffread(path1);
nFrame = length(dat1);
maxVal = 2^dat1(1).bits-1;
if iscell(dat1(1).data)
    [Nx,Ny] = size(dat1(1).data{1});
else
    [Nx,Ny] = size(dat1(1).data);
end

sigch1 = zeros(Nx,Ny,nFrame);
for ii=1:nFrame
    tmp = dat1(ii);
    if iscell(tmp.data)
        datEle = double(tmp.data{1})/maxVal;
    else
        datEle = double(tmp.data)/maxVal;
    end
    if ii==1
        sigch1Mean = datEle;
    else
        sigch1Mean = sigch1Mean + datEle;
    end
    sigch1(:,:,ii) = datEle;
end
sigch1Mean = sigch1Mean/nFrame;

if crop
    dat = sigch1(500:end,500:end,:);
    datRef = sigch1Mean(500:end,500:end);
else
    dat = sigch1;
    datRef = sigch1Mean;
end
dats = reshape(dat,[],nFrame);

end