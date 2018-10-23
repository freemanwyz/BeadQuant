function [K0,occupyMap] = myAddText(K0, idx, ts, col, occupyMap)
% MYADDTEXT add text to image
% K0: input and output image
% idx: position to add the text in the image
% ts: digits
% d5: digit symbol template, from 0 to 9
% col: RGB color [a b c], all from 0 to 1

[Nx,Ny,~] = size(K0);
% occupyMap(sub2ind([Nx,Ny],idx(:,1),idx(:,2))) = 1;  % This should be done adding label

% digits
digitTemp = load('./cfg/digitTemp.mat');
d5 = digitTemp.res5/255*0.99;
[ndx,ndy,~] = size(d5);  % x is height and y is width

% if ts==3 || ts==15319
%     keyboard
% end

ts = num2str(ts);
nTs = length(ts);
tL = ndy*nTs;
tH = ndx;
tL2 = round(tL/2);

%% position
x0 = round(mean(idx(:,1)));
y0 = round(mean(idx(:,2)));
% x0 = idx(1,1);  % along height
% y0 = idx(1,2);  % along width
% ofst = 5;
ofst = round(max(abs(x0 - idx(:,1))));
bofst = 1;

% Eight candidate positions
xs = [-tH-ofst, -tH-ofst, -tH-ofst, 0,        0,    ofst,     ofst, ofst];
ys = [-tL-ofst, -tL2,     ofst,     -tL-ofst, ofst, -tL-ofst, -tL2, ofst];
occuNum = zeros(1,8);
% iirg = [2 4 5 7 1 3 6 8];
iirg = 1:8;
for ii=iirg
    x1 = x0 + xs(ii);
    y1 = y0 + ys(ii);
    xrg1 = x1:(x1+tH+bofst);
    yrg1 = y1:(y1+tL+bofst);
    % should within the image
    if max(xrg1)>=Nx || max(yrg1)>=Ny || min(xrg1)<1 || min(yrg1)<1
        occuNum(ii) = 1e9;
    else
        occ1 = occupyMap(xrg1,yrg1);
        occuNum(ii) = sum(occ1(:));
    end
end

% choose the least occupied valid position
[~,I]=min(occuNum);
x1 = x0 + xs(iirg(I));
y1 = y0 + ys(iirg(I));

% if y0+textLen+25 < Ny
%     y1 = y0 + 15;
% else
%     y1 = y0 - 15 - textLen;
% end
% if x0+textHt+25 < Nx
%     x1 = x0 + 15;
% else
%     x1 = x0 - 15 - textHt;
% end

%% add text with color
rgdx = x1:(x1+ndx-1);
rgdy = y1:(y1+ndy-1);
for kk=1:nTs
    idxx = str2double(ts(kk))+1;
    ly1 = K0(rgdx,rgdy,1);
    ly2 = K0(rgdx,rgdy,2);
    ly3 = K0(rgdx,rgdy,3);
    mk0 = d5(:,:,idxx);
    mk0Idx = mk0>0;
    ly1(mk0Idx) = mk0(mk0Idx)*col(1);
    ly2(mk0Idx) = mk0(mk0Idx)*col(2);
    ly3(mk0Idx) = mk0(mk0Idx)*col(3);
    K0(rgdx,rgdy,1) = ly1;
    K0(rgdx,rgdy,2) = ly2;
    K0(rgdx,rgdy,3) = ly3;
%     K0(rgdx,rgdy,1) = d5(:,:,idxx)*col(1);
%     K0(rgdx,rgdy,2) = d5(:,:,idxx)*col(2);
%     K0(rgdx,rgdy,3) = d5(:,:,idxx)*col(3);
    occupyMap(rgdx,rgdy) = occupyMap(rgdx,rgdy) + 1;
    xa = (min(rgdx)-5):(max(rgdx)+5);
    ya = (min(rgdy)-5):(max(rgdy)+5);
    occupyMap(xa,ya) = occupyMap(xa,ya) + 0.1;
    rgdy = rgdy + ndy;
end

end