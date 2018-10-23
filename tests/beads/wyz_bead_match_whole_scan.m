%% Bead signal from grace
p = getInitParamBead();
p.crop = 1;

myRg = 4:0.5:13;
myRef = beadTemplate(myRg);

% p.crop_rgx = 1851:2000;
% p.crop_rgy = 1751:1900;

% p.crop_rgx = 419:647;
% p.crop_rgy = 48:200;

% p.crop_rgx = 80:230;
% p.crop_rgy = 1900:2200;

% p.crop_rgx = 1501:1700;
% p.crop_rgy = 1301:1600;

% p.crop_rgx = 651:800;
% p.crop_rgy = 1051:1150;

% p.crop_rgx = 3178:3328;
% p.crop_rgy = 1:1500;

% load data -----
fin = [p.tp filesep 'dump' filesep 'beadAvg.tif'];
dat = double(imread(fin));

[NxA,NyA] = size(dat);
gX = 1:100:NxA;
gY = 1:100:NyA;
for ii=gX
    for jj=gY
        p.crop_rgx = ii:(ii+125);
        p.crop_rgx = p.crop_rgx( p.crop_rgx<=NxA );
        p.crop_rgy = jj:(jj+125);
        p.crop_rgy = p.crop_rgy( p.crop_rgy<=NyA );
        
        % To isolated parts ----
        if p.crop
            dat0 = dat(p.crop_rgx,p.crop_rgy);
        else
            dat0 = dat;
        end
        % dat0 = log(dat0+1);
        % xMax = log(255+1);
        % xMax = log(255+sqrt(255^2+1));
        
%         xMax = log(1+sqrt(1^2+1));
        xMax = 1;
        
        dat0 = dat0/255;
        
%         dat0a = log(dat0+sqrt(dat0.^2+1));
        dat0a = sqrt(dat0);
        
        dat0a = dat0a/xMax;
        [Nx,Ny] = size(dat0);
        
        % detecting ---
        resBead = beadMatchPairSeqNoThr( dat0a, myRef, 0 );
        % resBead = beadMatchWhole( dat0, myRef, 0 );
        
        % plot ---
        resBorder = dat0*0;
        neibVec = [0, -1, 1, -Nx, Nx];
        
        fprintf('Dump ===== \n');
        for tt=1:length(resBead)
            idx = resBead{tt};
            idx = sub2ind(size(dat0),idx(:,1),idx(:,2));
            idxk = reshape(bsxfun(@plus,idx,neibVec),[],1);
            idxk = idxk(idxk>0 & idxk<=(Nx*Ny));
            idxc = setdiff(idxk,idx);
            resBorder(idxc) = 0.75;
        end
        
        K2 = cat(3,resBorder,dat0a,dat0a*0);
%         K2 = cat(3,resBorder,sqrt(dat0),dat0*0);
        fname = [p.tp 'dump' filesep '20150627_sqrt' filesep num2str(ii) '_' num2str(jj) '.tif'];
        imwrite(K2,fname);
%         figure;imshow(K2);
    end
end








