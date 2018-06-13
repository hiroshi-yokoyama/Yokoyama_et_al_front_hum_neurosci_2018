%    x : input data: 2D or 3D matrix [ch x time] / [freqs x time] / [ch x time x trial] / [freqs x time x trial]
% indx : baseline index of time [1 x sample]
function z_value = z_score (x, indx)
    if length(size(x))==2
        MEAN = mean(x(:,indx),2);
        SD   = std(x(:,indx),0,2);        
    elseif length(size(x))==3
        MEAN = mean(x(:,indx,:),2);
        SD   = std(x(:,indx,:),0,2);
    elseif length(size(x))==4
        MEAN = mean(x(:,indx,:,:),2);
        SD   = std(x(:,indx,:,:),0,2);
    elseif length(size(x))==5
        MEAN = mean(x(:,indx,:,:,:),2);
        SD   = std(x(:,indx,:,:,:),0,2);
    end
    
    z_value = bsxfun(@rdivide, bsxfun(@minus, x, MEAN), SD);
end



















