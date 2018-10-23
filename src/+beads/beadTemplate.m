function res = beadTemplate( radRg, Nx )
% BEADTEMPLATE Build pairs around the border of a cicle

% Nx = 41;
ofstCenter = (Nx-1)/2+1;
neibVec = [-1,1,-Nx,Nx,-Nx-1,-Nx+1,Nx-1,Nx+1];
% neibVec4 = [-1,1,-Nx,Nx];

res = [];
nRad = length(radRg);
for ii=1:nRad
    tmp = zeros(Nx,Nx);
    [ax,ay] = find(tmp>-100);
    C0 = radRg(ii)-0.01;  % to make the circle look better
    % for each one in aOut0, search suitable aOut1 and aOut2 -----
    aDist = sqrt((ax-ofstCenter).^2 + (ay-ofstCenter).^2);
    
    % full circle plus one
%     aOutAll = aDist <= C0;
    aOutAll = aDist <= C0 + 1;
    tmpAll = zeros(Nx,Nx);
    tmpAll(aOutAll) = 1;
    [ax,ay] = find(tmpAll > 0);
    
    % find a ring
    aOut0 = aDist>C0-1 & aDist<=C0;
    tmp(aOut0) = 1;

    refPix = reshape(find(aOut0),[],1);    
    nRefPix = length(refPix);
    myRef = zeros(nRefPix,2+1);
    myRefX = zeros(nRefPix,2+1);
    myRefY = zeros(nRefPix,2+1);
    
    for jj=1:nRefPix
        pixMe = refPix(jj);
        [pixMeX,pixMeY] = ind2sub(size(tmp),pixMe);
        angleOut = angle(complex(pixMeX-ofstCenter, pixMeY-ofstCenter))*180/pi;
        angleIn = angle(complex(ofstCenter-pixMeX, ofstCenter-pixMeY))*180/pi;
        pixNeib = pixMe + neibVec;
        [pixNeibX,pixNeibY] = ind2sub(size(tmp),pixNeib);
        angleNeib = angle(complex(pixNeibX-pixMeX,pixNeibY-pixMeY))*180/pi;
        angleDifOut = abs(angleNeib - angleOut);
        angleDifIn = abs(angleNeib - angleIn);
        [~,outI] = min(angleDifOut);
        [~,inI] = min(angleDifIn);
        myRef(jj,:) = [pixMe,pixNeib(inI),pixNeib(outI)];
        myRefX(jj,:) = [pixMeX,pixNeibX(inI),pixNeibX(outI)];
        myRefY(jj,:) = [pixMeY,pixNeibY(inI),pixNeibY(outI)];
    end
%     K0 = tmp*0;
%     K0(myRef(:,1)) = 1;
%     K0(myRef(:,2)) = 0.5;
%     K0(myRef(:,3)) = 0.25;
    
    res(ii).mask = myRef;
    res(ii).maskX = myRefX;
    res(ii).maskY = myRefY;
    res(ii).radius = radRg(ii);
    res(ii).maskCir = [ax,ay];
end





