function bounds = mask2bounds(mask)
%mask2bounds: takes a 1d logical vector and returns the boundaries of all
%strings of trues.
%   mask: 1d logical array
%   bounds: nx2 array of boundaries for n regions

if ~islogical(mask)
    %error('mask must be a logical vector')
    mask=logical(mask);
end
if size(mask,1) == 1
    if size(mask,2) > 1
        mask = mask';
    end
elseif size(mask,2) > 1
    error('mask must be 1 dimensional')
end

on = find(diff(mask)==1) + 1;
off = find(diff(mask)==-1);

if mask(1) == true
    on = [1; on];
end
if mask(end) == true
    off = [off; length(mask)];
end

bounds = [on off];

end

