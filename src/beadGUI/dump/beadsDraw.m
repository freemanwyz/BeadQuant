function beadsDraw(resBead,bgMean,sig,goodBeadIdx,controlBeadIdx,outPath,maxVal)
% BEADDRAW Draw beads on bgMean (background average)
% resBead: coordinates for all beads
% goodBeadIdx: index of beads to plot

[Nx,Ny] = size(bgMean);
xI = goodBeadIdx;
neibVec = [0, -1, 1, -Nx, Nx];

%% plot all selected beads on the averaged background
t0 = zeros(Nx,Ny);      % selected beads
for ii=1:length(xI)
    idx = resBead{xI(ii)};
    idx = sub2ind(size(t0),idx(:,1),idx(:,2));
    idxk = bsxfun(@plus,idx,neibVec);
    idxa = ismember(idxk,idx);
    idxSel = sum(idxa,2)<5;
    idxc = idx(idxSel);
    t0(idxc) = 1;
end
t1 = zeros(Nx,Ny);      % control beads
for ii=1:length(controlBeadIdx)
    idx = resBead{controlBeadIdx(ii)};
    idx = sub2ind(size(t0),idx(:,1),idx(:,2));
    idxk = bsxfun(@plus,idx,neibVec);
    idxa = ismember(idxk,idx);
    idxSel = sum(idxa,2)<5;
    idxc = idx(idxSel);
    t1(idxc) = 1;
end
tr = bgMean*0;
tg = bgMean;
tb = tr;
tr(t1>0 | t0>0) = 1;
tg(t1>0 | t0>0) = 1;
tb(t0>0) = 1;
K0 = cat(3,tr,tg,tb);
% imshow(K0);
% f0 = figure('visible','off');
% imshow(K0);

%% label beads
% avoid overlapping with selected or control beads
K_occupy = zeros(Nx,Ny);
for ii=1:length(xI)
    idx = resBead{xI(ii)};
    K_occupy(sub2ind([Nx,Ny],idx(:,1),idx(:,2))) = 5;
end
for ii=1:length(controlBeadIdx)
    idx = resBead{controlBeadIdx(ii)};
    K_occupy(sub2ind([Nx,Ny],idx(:,1),idx(:,2))) = 5;
end

% add label
for ii=1:length(xI)
    idx = resBead{xI(ii)}; 
    [K0, K_occupy] = myAddText(K0, idx, xI(ii), [1 1 1], K_occupy);    
%     K0 = insertText(K0,a0,num2str(xI(ii)),'TextColor',...
%         [0.99 0.99 0.99],'FontSize',50,'BoxOpacity',0.75,'BoxColor','black');
%     text(idx(1,2)+15,idx(1,1)+15,num2str(xI(ii)),'Color',[0.99 0.99 0.99],'FontSize',10);
end
for ii=1:length(controlBeadIdx)
    idx = resBead{controlBeadIdx(ii)};
    [K0,K_occupy] = myAddText(K0, idx, controlBeadIdx(ii), [1 1 0], K_occupy);    
%     K0 = insertText(K0,a0,num2str(controlBeadIdx(ii)),...
%         'TextColor',[0.99 0.99 0],'FontSize',50,'BoxOpacity',0.75,'BoxColor','black');
%     text(idx(1,2)+15,idx(1,1)+15,num2str(controlBeadIdx(ii)),'Color',[1 1 0],'FontSize',10);
end

% display and save
imshow(K0);
% set(f0,'PaperPositionMode','auto')
% set(gca,'Position',[0.05 0.05 0.9 0.9])
h = msgbox('Saving image...');
% print(f0,'-dpng','-r2500',[outPath,filesep,'bead_selected.png'])
f00 = [outPath,filesep,'bead_selected.png'];
imwrite(double(K0),f00);
if isvalid(h)
    close(h);
end
% close(f0);

%% plot selected beads on signal channel one by one
h = msgbox('Saving figures...');
nTps = size(sig,3);
for ii=1:length(xI)
    f1 = [outPath,filesep,'bead_',num2str(xI(ii)),'_signal_channel.tif'];
    pix0 = resBead{xI(ii)};
    Row0 = round((min(pix0(:,1)) + max(pix0(:,1)))/2);
    Col0 = round((min(pix0(:,2)) + max(pix0(:,2)))/2);
    xrg = (Row0-100):(Row0+100);
    yrg = (Col0-100):(Col0+100);
    xrg = xrg(xrg>0 & xrg<=Nx);
    yrg = yrg(yrg>0 & yrg<=Ny);
    t0Ele = t0(xrg,yrg);
    for tt=1:nTps        
        if tt==1
            K1 = cat(3,t0Ele*0.5,t0Ele*0,double(sig(xrg,yrg,tt))/maxVal);
            imwrite(double(K1),f1);
        else
            K1 = cat(3,t0Ele*0.5,t0Ele*0,double(sig(xrg,yrg,tt))/maxVal);
            imwrite(double(K1),f1,'WriteMode','append');
        end
    end
end

f1 = [outPath,filesep,'bead_selected_signal_channel.tif'];
if exist(f1, 'file')==2
  delete(f1);
  ss = [f1 ' already exist\n'];
  warning(ss);
end
for tt=1:nTps
    K1 = cat(3,t0,t0*0,double(sig(:,:,tt))/maxVal);
    if tt==1
        imwrite(double(K1),f1);
    else
        imwrite(double(K1),f1,'WriteMode','append');
    end
    pause(0.25);
end

%% plot all beads, highlight selected beads
resBorder = zeros(Nx,Ny);
for ii=1:length(resBead)
    if mod(ii,1000)==0
        fprintf('ii is %d\n',ii);
    end
    idx = resBead{ii};
    idx = sub2ind(size(resBorder),idx(:,1),idx(:,2));
    idxk = bsxfun(@plus,idx,neibVec);
    idxa = ismember(idxk,idx);
    idxSel = sum(idxa,2)<5;
    idxc = idx(idxSel);
    resBorder(idxc) = 0.5;
end
K2 = cat(3,resBorder.*(1-t0),bgMean.*(1-t0),t0);
fname = 'res_bead_highlight';
imwrite(double(K2),[outPath,filesep,fname,'.tif']);

% plot all beads, highlight selected beads and control beads
resBorder = zeros(Nx,Ny);
for ii=1:length(resBead)
    if mod(ii,1000)==0
        fprintf('ii is %d\n',ii);
    end
    idx = resBead{ii};
    idx = sub2ind(size(resBorder),idx(:,1),idx(:,2));
    idxk = bsxfun(@plus,idx,neibVec);
    idxa = ismember(idxk,idx);
    idxSel = sum(idxa,2)<5;
    idxc = idx(idxSel);
    resBorder(idxc) = 0.5;
end
ta = resBorder.*(1-t0) + t1;
ta(ta>1) = 1;
tb = bgMean.*(1-t0);
tb(t1>0) = 1;
K3 = cat(3,ta,tb,t0);
fname = 'res_bead_highlight_with_control';
imwrite(double(K3),[outPath,filesep,fname,'.tif']);

if isvalid(h)
    close(h);
end

end


