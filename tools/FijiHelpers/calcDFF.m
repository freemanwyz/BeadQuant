function DFF_mat = calcDFF(F_mat,timeDim,MODEFLAG,baseline)
%This function calculates a DFF matrix with time along columns. User must
%input the time dimension of input matrix (timeDim) and the baseline period
%over which to average F0 unless MODEFLAG is set to 1.
%
%If MODEFLAG is set to 1, values within 1 SD of the mode are averaged to
%calculate F0. This is useful in cases where heightened fluorescence is
%sparse as is typical for, e.g. spontaneous astrocyte calcium data.

if timeDim~=1 %if time dose not lie along columns, take matrix transpose
    F_mat = F_mat';
end

if ~MODEFLAG
    %Calc DFF based on declared baseline SD
    F0_vect = mean(F_mat(1:baseline,:),1);
    F0_mat = repmat(F0_vect,[size(F_mat,1),1]);
    DFF_mat = (F_mat - F0_mat)./F0_mat;
    
elseif MODEFLAG
    %Get upper and lower limits to average for DFF
    F_mode = mode(F_mat,1);
    F_SD = std(F_mat,0,1);
    FU_vect = F_mode + F_SD;%vector of values one SD above each trace mode
    FL_vect = F_mode - F_SD;%vector of values one SD below each trace mode
    
    for iROI = 1:size(F_mat,2)
        F = F_mat(:,iROI);
        F_good = F(F<=FU_vect(iROI) & F>=FL_vect(iROI));%F values within range
        F0 = mean(F_good);
        DFF_mat(:,iROI) = (F - F0)/F0;
    end
    
end