function [fmat,ROIsize_vect] = extractFmat
%This function extracts the fluorescence time series data from matrices
%produced by fiji multi-measure tool when applied to ROIs. It returns these
%as column vectors packaged into a matrix as well as a vector containing
%the size of each ROI
fijiDat = paste;
fmat = fijiDat(:,3:4:end);
ROIsize_vect = fijiDat(1,2:4:end);