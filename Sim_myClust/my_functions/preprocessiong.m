% Preprocessing function for EEG signals
%  [data]=preprocessiong(raw_data,base_index,filtering,zscore)
%  raw_data : 補正前の生データ ch x times
%  base_index : row or colmn vector of baseline index. 
%               e.g.) base_index=1:512
%  filtering  : cutoff frequency parameter
%              If set as "filtering = [1 50]./(1024/2)", Butterworth filter with cutoff frequency between 1 to 50 Hz is applied. 
%              If set as "[]", filtering is skiped for preprocessing         
%  zscore     : z-score flag
%               If set as "on", z-score is applied with mean and SD valuse in baseline periods
function [outdata]=preprocessiong(raw_data,base_index,filtering,zscore)
    raw_data=double(raw_data);
    data=raw_data;
    
    base_data=data(:,base_index);
    mean_base=mean(base_data,2);
    data=bsxfun(@minus, data, mean_base);
    
    if ~isempty(filtering)
        [b,a]=butter(3,filtering);
        data=filtfilt(b,a,data.').';        
    end
    
    base_data=data(:,base_index);
    mean_base=mean(base_data,2);
    data=bsxfun(@minus, data, mean_base);
    
    if strcmp(zscore,'on')
        std_base=std(base_data,0,2);
        data=bsxfun(@rdivide, data, std_base);
    end
    
    outdata=data;
    
end