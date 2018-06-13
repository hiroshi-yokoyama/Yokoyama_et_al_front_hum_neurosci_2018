% tf     : time frequency spectrum [freqs x time x trial x ch]
% wpli   : estimated value of weighted phase lag indecies
% mode   : mode=0; normal wpli, mode=1; debiased wpli
function wpli = cal_wpli(tf, mode)
    if nargin<2
        mode=0; % normal wpli
    end
    [Nfrq, Nt, Ntri, Nch]=size(tf);
    
    wpli=zeros(Nfrq, Nt, Nch, Nch);
    
    for ch1=1:Nch
        for ch2=1:Nch
            X1=tf(:,:,:,ch1);
            X2=tf(:,:,:,ch2);
            
            X12 = X1.*conj(X2);
            imX12    = imag(X12);

            outsum   = nansum(imX12,3);      % compute the sum
            outsumW  = nansum(abs(imX12),3); % normalization of the WPLI

            outssq   = nansum(imX12.^2,3);
            if mode ==0
                tmp_wpli     = (outsum)./outsumW;
            elseif mode==1 % debiased method
                tmp_wpli     = (outsum.^2 - outssq)./(outsumW.^2 - outssq);%
            end
            
            wpli(:, :, ch1, ch2)=abs(tmp_wpli);
            
            disp(['ch', num2str(ch1), ' vs ', 'ch', num2str(ch2), ' : OK'])
        end
    end
end