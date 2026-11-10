% Used for: nitracline depth from a nitrate profile (used by zno3_3threshold_run.m)
% Original: Code\CCHDO\get_nitracline.m (data paths updated to I:\Prj; computation unchanged)
function nitracline_depth = get_nitracline(no3_profile, depth_vec, threshold)
%GET_NITRACLINE_DEPTH  Compute nitracline depth (NO3 = threshold) for
%multiple profiles using linear interpolation.
%
%   nitracline_depth = GET_NITRACLINE_DEPTH(no3_profile, depth_vec, threshold)
%
%   Inputs:
%       no3_profile : [Ndepth × Nprofiles] matrix
%       depth_vec   : [Ndepth × 1] depth vector (positive downward)
%       threshold   : scalar threshold (e.g., 3 µmol/kg)
%
%   Output:
%       nitracline_depth : [Nprofiles × 1] vector of nitracline depths

    [Ndepth, Nprof] = size(no3_profile);
    depth_vec = depth_vec(:);  % ensure column vector

    nitracline_depth = NaN(Nprof,1);  % output

    for p = 1:Nprof
        x = no3_profile(:,p);

        % First depth index where NO3 >= threshold
        idx = find(x >= threshold, 1, 'first');

        if isempty(idx)
            % Never reaches threshold → leave NaN
            continue;
        end

        if idx == 1
            nitracline_depth(p) = depth_vec(1);
            continue;
        end

        % Interpolate between the two levels
        x1 = x(idx-1);  x2 = x(idx);
        y1 = depth_vec(idx-1);  y2 = depth_vec(idx);

        if isnan( x1 )
            nitracline_depth(p) = nan;
        else
            nitracline_depth(p) = interp1([x1 x2], [y1 y2], threshold);
        end
    end
end
