function data = remove_speckle(data)
% Function to remove speckle. 
% First quick version.
% RG 29.10.2020

mask = any(~isnan( data.spec),3);

mask1 = repmat( mask , [1 1 max(data.DoppLen)]);

for tt = 1:data.totsamp
    for rr = data.range_offsets(end):data.n_levels
%    for rr = 1:data.n_levels
    
        if mask(tt,rr) == 0; continue; end
        
        % deal with edges
        if tt == 1
            ts = tt;
            te = tt+1;
        elseif tt == data.totsamp
            ts = tt-1;
            te = tt;
        else
            ts = tt-1;
            te = tt+1;
        end
        
        if rr == 1
            rs = rr;
            re = rr+1;
        elseif rr == data.n_levels
            rs = rr-1;
            re = rr;
        else
            rs = rr-1;
            re = rr+1;
        end
        
        if sum(mask(ts:te,rs:re)) <= 2;
            mask1(tt,rr,:) = 0;
        end
        
        %if sum(mask(tt-1:tt+1,rr-1:rr+1)) <= 2;
        %    mask1(tt,rr,:) = 0;
        %end
        
    end % rr
end % tt    


mask1 = ~mask1;    
data.spec(mask1) = NaN;    

if isfield(data,'spec_hv')
    data.spec_hv(mask1) = NaN;    
end    


if isfield(data,'spec_covRe')
    data.spec_covRe(mask1) = NaN;
    data.spec_covIm(mask1) = NaN;
end
