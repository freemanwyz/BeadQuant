function [accumMatrix, gradientImg] = chaccum(A,radiusRange,mask)
%CHACCUM Compute 2D accumulator array using Circular Hough Transform

maxNumElemNHoodMat = 1e7; % Maximum number of elements in neighborhood matrix xc allowed, before memory chunking kicks in.

%% Check if the image is flat
flat = all(A(:) == A(1));
if (flat)
    accumMatrix = zeros(size(A,1),size(A,2));
    gradientImg = zeros(size(A,1),size(A,2));
    return;
end

%% Get the input image in the correct format
A = getGrayImage(A);

%% Calculate gradient
% with mask
[Gx,Gy,gradientImg] = imgradient(A);
gradientImg = gradientImg.*mask;

%% Get edge pixels
[Ex, Ey] = getEdgePixels(gradientImg, -1);
idxE = sub2ind(size(gradientImg), Ey, Ex);

%% Identify different radii for votes
if (length(radiusRange) > 1)
    radiusRange = radiusRange(1):0.5:radiusRange(2);
end
RR = radiusRange;

%% Compute the weights for votes for different radii
if (length(radiusRange) > 1)
    lnR = log(radiusRange);
    phi = ((lnR - lnR(1))/(lnR(end) - lnR(1))*2*pi) - pi; % Modified form of Log-coding from Eqn. 8 in [3]
else
    phi = 0;
end
Opca = exp(sqrt(-1)*phi);
w0 = Opca./(2*pi*radiusRange);

%% Computing the accumulator array
xcStep = floor(maxNumElemNHoodMat/length(RR));
lenE = length(Ex);
[M, N] = size(A);
accumMatrix = zeros(M,N);

for i = 1:xcStep:lenE
    Ex_chunk = Ex(i:min(i+xcStep-1,lenE));
    Ey_chunk = Ey(i:min(i+xcStep-1,lenE));
    idxE_chunk = idxE(i:min(i+xcStep-1,lenE));
    
    xc = bsxfun(@plus, Ex_chunk, bsxfun(@times, -RR, Gx(idxE_chunk)./gradientImg(idxE_chunk))); % Eqns. 10.3 & 10.4 from Machine Vision by E. R. Davies
    yc = bsxfun(@plus, Ey_chunk, bsxfun(@times, -RR, Gy(idxE_chunk)./gradientImg(idxE_chunk)));
    
    xc = round(xc);
    yc = round(yc);
    
    w = repmat(w0, size(xc, 1), 1);
    
    % weighted circular Hough transform
    val0 = gradientImg(sub2ind(size(gradientImg),Ey_chunk,Ex_chunk));
    val1 = repmat(val0,1,size(w,2));    
    w = w.*val1;
    
    %% Determine which edge pixel votes are within the image domain
    % Record which candidate center positions are inside the image rectangle.
    [M, N] = size(A);
    inside = (xc >= 1) & (xc <= N) & (yc >= 1) & (yc < M);
    
    % Keep rows that have at least one candidate position inside the domain.
    rows_to_keep = any(inside, 2);
    xc = xc(rows_to_keep,:);
    yc = yc(rows_to_keep,:);
    w = w(rows_to_keep,:);
    inside = inside(rows_to_keep,:);
    
    %% Accumulate the votes in the parameter plane
    xc = xc(inside); yc = yc(inside);
    win = w(inside);
    accumMatrix = accumMatrix + accumarray([yc(:), xc(:)], win, [M, N]);
    clear xc yc w; % These are cleared to create memory space for the next loop. Otherwise out-of-memory at xc = bsxfun... in the next loop.
end

end

function [Gx, Gy, gradientImg] = imgradient(I)
hy = -fspecial('sobel');
hx = hy';
Gx = imfilter(I, hx, 'replicate','conv');
Gy = imfilter(I, hy, 'replicate','conv');
if nargout > 2
    gradientImg = hypot(Gx, Gy);
end
end

function [Ex, Ey] = getEdgePixels(gradientImg, edgeThresh)
Gmax = max(gradientImg(:));
if (isempty(edgeThresh))
    edgeThresh = graythresh(gradientImg/Gmax); % Default EdgeThreshold
end
t = Gmax * edgeThresh;
[Ey, Ex] = find(gradientImg > t);
end

function A = getGrayImage(A)
N = ndims(A);
if (N == 3) % RGB Image
    A = rgb2gray(A);
    if (isinteger(A))
        A = im2single(A); % If A is an integer, cast it to floating-point
    end    
elseif (N == 2)
    if (islogical(A)) % Binary image
        filtStd = 1.5;
        filtSize = ceil(filtStd*3);
        filtSize = filtSize + ceil(rem(filtSize,2)); % filtSize = Smallest odd integer greater than filtStd*3
        gaussFilt = fspecial('gaussian',[filtSize filtSize],filtStd);
        A = imfilter(im2single(A),gaussFilt,'replicate');
    elseif (isinteger(A))
        A = im2single(A); % If A is an integer, cast it to floating-point
    end    
else
    iptassert(false,'images:imfindcircles:invalidInputImage'); % This should never happen here.
end

end

