function [centers, metric] = chcenters(accumMatrix,suppThreshold)
%CHCENTERS Find circle center locations from the Circular Hough Transform accumulator array 

sigma = [];
medFiltSize = 3; % Size of the median filter
centers = [];
metric = [];
%% Use the magnitude - Accumulator array can be complex. 
accumMatrix = abs(accumMatrix);

%% Check if the accumulator array is flat
flat = all(accumMatrix(:) == accumMatrix(1));
if (flat)
    return;
end
%% Filter the accumulator array
if (~isempty(sigma))
    accumMatrix = gaussianFilter(accumMatrix, sigma);
end

%% Pre-process the accumulator array
if (min(size(accumMatrix)) > medFiltSize)
    Hd = medfilt2(accumMatrix, [medFiltSize medFiltSize]); % Apply median filtering only if the image is big enough.
else
    Hd = accumMatrix;
end
suppThreshold = max(suppThreshold - eps(suppThreshold), 0);
Hd = imhmax(Hd, suppThreshold);
bw = imregionalmax(Hd);
s = regionprops(bw,accumMatrix,'weightedcentroid'); % Weighted centroids of detected peaks.

%% Sort the centers based on their accumulator array value
if (~isempty(s))
    centers = reshape(cell2mat(struct2cell(s)),2,length(s))';
    % Remove centers which are NaN.
    [rNaN, ~] = find(isnan(centers));
    centers(rNaN,:) = [];
    
    if(~isempty(centers))
        metric = Hd(sub2ind(size(Hd),round(centers(:,2)),round(centers(:,1))));
        % Sort the centers in descending order of metric
        [metric,sortIdx] = sort(metric,1,'descend');
        centers = centers(sortIdx,:);
    end
end

end

function accumMatrix = gaussianFilter(accumMatrix, sigma)
    filtSize = ceil(sigma*3);
    filtSize = min(filtSize + ceil(rem(filtSize,2)), min(size(accumMatrix))); % filtSize = Smallest odd integer greater than sigma*3
    gaussFilt = fspecial('gaussian',[filtSize filtSize],sigma);
    accumMatrix = imfilter(accumMatrix, gaussFilt, 'same');
end




