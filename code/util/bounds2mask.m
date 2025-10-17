function mask = bounds2mask(bounds, mask_len, num_pad, num_pad_end)
%BOUNDS2MASK takes an n x 2 matrix representing the bounds of n regions in
%a 1d array and returns a logical mask of length mask_len representing
%those regions, with each region optionally padded on both sides by num_pad
%elements, or if 4 args, padded before the region by num_pad and after by
%num_pad_end
% bounds2mask(bounds, mask_len, num_pad, num_pad_end)
if nargin < 3
    num_pad = 0;
end
if nargin < 4
    num_pad_end = num_pad;
end

mask = false(mask_len,1);
if ~isempty(bounds)
    for n = 1:size(bounds,1)
        curr_bounds = bounds(n,:);
        curr_bounds(1) = curr_bounds(1) - num_pad;
        curr_bounds(2) = curr_bounds(2) + num_pad_end;
        curr_bounds(curr_bounds<1) = 1;
        curr_bounds(curr_bounds>mask_len) = mask_len;
        if(curr_bounds(1)<=curr_bounds(2))
            mask(curr_bounds(1):curr_bounds(2)) = true;
        end
    end
end
end

