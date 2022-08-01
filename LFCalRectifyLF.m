% LFCalRectifyLF - rectify a light field using a calibrated camera model, called as part of LFUtilDecodeLytroFolder
%
% Usage:
%     [LF, RectOptions] = LFCalRectifyLF( LF, CalInfo, RectOptions )
%     [LF, RectOptions] = LFCalRectifyLF( LF, CalInfo )
%
% This function is called by LFUtilDecodeLytroFolder to rectify a light field. It follows the
% rectification procedure described in:
%
% D. G. Dansereau, O. Pizarro, and S. B. Williams, "Decoding, calibration and rectification for
% lenslet-based plenoptic cameras," in Computer Vision and Pattern Recognition (CVPR), IEEE
% Conference on. IEEE, Jun 2013.
%
% Minor differences from the paper: light field indices [i,j,k,l] are 1-based in this
% implementation, and not 0-based as described in the paper.
%
% Note that a calibration should only be applied to a light field decoded using the same lenslet
% grid model. This is because the lenslet grid model forms an implicit part of the calibration,
% and changing it will invalidate the calibration.
%
% LFCalDispRectIntrinsics is useful for visualizing the sampling pattern associated with a requested
% intrinsic matrix.
%
% Inputs:
%
%         LF : The light field to rectify; should be floating point
%
%     CalInfo struct contains a calibrated camera model, see LFUtilCalLensletCam:
%            .EstCamIntrinsicsH : 5x5 homogeneous matrix describing the lenslet camera intrinsics
%            .EstCamDistortionV : Estimated distortion parameters
%
%    [optional] RectOptions struct (all fields are optional) :
%        .NInverse_Distortion_Iters : Number of iterations in inverse distortion estimation
%                        .Precision : 'single' or 'double'
%               .RectCamIntrinsicsH : Requests a specific set of intrinsics for the rectified light
%                                     field. By default the rectified intrinsic matrix is
%                                     automatically constructed from the calibrated intrinsic
%                                     matrix, but this process can in some instances yield poor
%                                     results: excessive black space at the edges of the light field
%                                     sample space, or excessive loss of scene content off the edges
%                                     of the space. This parameters allows you to fine-tune the
%                                     desired rectified intrinsic matrix.
%
% Outputs :
%
%     LF : The rectified light field
%     RectOptions : The rectification options as applied, including any default values employed.
%
% User guide: <a href="matlab:which LFToolbox.pdf; open('LFToolbox.pdf')">LFToolbox.pdf</a>
% See also: LFCalDispRectIntrinsics, LFUtilCalLensletCam, LFUtilDecodeLytroFolder, LFUtilProcessWhiteImages

% Copyright (c) 2013-2020 Donald G. Dansereau

function [LF, RectOptions] = LFCalRectifyLF(LF, CalInfo, RectOptions, CachedIdx)

    %---Defaults---
    RectOptions = LFDefaultField('RectOptions', 'MaxUBlkSize', 32);
    RectOptions = LFDefaultField('RectOptions', 'Precision', 'single');
    LFSize = size(LF);

    %---Build interpolation indices---
    fprintf('Generating interpolation indices...\n');
    LF = cast(LF, RectOptions.Precision);

    fprintf('Interpolating...');
    tStart = cputime;
    LFOut = LF;

    if exist('CachedIdx', 'var') && ~isempty(CachedIdx)
        % Use the cached index values

        % Interpolate the entire light field in one go
        LFOut = InterpolateColours(LF, LFOut, CachedIdx, LFSize(1:4), 1:LFSize(4));

    else
        % Default option, recompute InterpIdx
        UBlkSize = RectOptions.MaxUBlkSize;

        %---chop up the LF along u---

        for (UStart = 1:UBlkSize:LFSize(4))
            UStop = UStart + UBlkSize - 1;
            UStop = min(UStop, LFSize(4));
            DestSize = [LFSize(1:3), length(UStart:UStop)];
            
            % InterpIdx initially holds the index of the desired ray, and is evolved through the application
            % of the inverse distortion model to eventually hold the continuous-domain index of the undistorted
            % ray, and passed to the interpolation step.
            [InterpIdx, RectOptions] = LFCalComputeIdx(LFSize, [arrayfun(@(x)1:x, LFSize(1:3), 'UniformOutput', false), UStart:UStop], ...
                CalInfo, RectOptions);
            LFOut = InterpolateColours(LF, LFOut, InterpIdx, DestSize, UStart:UStop);

            fprintf('.')
        end

    end

    LF = LFOut;
    clear LFOut;

    %---Clip interpolation result, which sometimes rings slightly out of range---
    LF(isnan(LF)) = 0;
    LF = max(0, min(1, LF));

    fprintf('\nDone. Interpolation took %.2f CPU seconds.\n', cputime - tStart);
end

function LFOut = InterpolateColours(LF, LFOut, InterpIdx, DestSize, USlice)
    LFSize = size(LF);
    NChans = LFSize(5);

    %{
        InterpIdx = min(LFSize(1:4)' - 1, max(1, InterpIdx));

        for (ColChan = 1:NChans)
            S0 = floor(InterpIdx(2,:));
            T0 = floor(InterpIdx(1,:));
            U0 = floor(InterpIdx(4,:));
            V0 = floor(InterpIdx(3,:));

            S1 = S0 + 1;
            T1 = T0 + 1;
            U1 = U0 + 1;
            V1 = V0 + 1;

            Sd = InterpIdx(2,:) - S0;
            Td = InterpIdx(1,:) - T0;
            Ud = InterpIdx(4,:) - U0;
            Vd = InterpIdx(3,:) - V0;

            GetLF = @(S, T, U, V, ColChan) LF(sub2ind(LFSize, S, T, U, V, repmat(ColChan, size(S))));

            c000 = GetLF(S0, T0, U0, V0, ColChan) .* (1 - Sd) + GetLF(S1, T0, U0, V0, ColChan) .* Sd;
            c001 = GetLF(S0, T0, U0, V1, ColChan) .* (1 - Sd) + GetLF(S1, T0, U0, V1, ColChan) .* Sd;
            c010 = GetLF(S0, T0, U1, V0, ColChan) .* (1 - Sd) + GetLF(S1, T0, U1, V0, ColChan) .* Sd;
            c011 = GetLF(S0, T0, U1, V1, ColChan) .* (1 - Sd) + GetLF(S1, T0, U1, V1, ColChan) .* Sd;
            c100 = GetLF(S0, T1, U0, V0, ColChan) .* (1 - Sd) + GetLF(S1, T1, U0, V0, ColChan) .* Sd;
            c101 = GetLF(S0, T1, U0, V1, ColChan) .* (1 - Sd) + GetLF(S1, T1, U0, V1, ColChan) .* Sd;
            c110 = GetLF(S0, T1, U1, V0, ColChan) .* (1 - Sd) + GetLF(S1, T1, U1, V0, ColChan) .* Sd;
            c111 = GetLF(S0, T1, U1, V1, ColChan) .* (1 - Sd) + GetLF(S1, T1, U1, V1, ColChan) .* Sd;

            c00 = c000 .* (1 - Td) + c100 .* Td;
            c01 = c001 .* (1 - Td) + c101 .* Td;
            c10 = c010 .* (1 - Td) + c110 .* Td;
            c11 = c011 .* (1 - Td) + c111 .* Td;

            c0 = c00 .* (1 - Ud) + c10 .* Ud;
            c1 = c01 .* (1 - Ud) + c11 .* Ud;

            c = c0 .* (1 - Vd) + c1 .* Vd;

            LFOut(:, :, :, :, ColChan) = reshape(c, LFSize(1:4));
        end
    %}

    % todo[optimization]: use a weighted interpolation scheme to exploit the weight channel

    % Using custom implentation of GriddedInterpolant may be faster. Custom interp1d exists on file exchange that
    % can be faster than interp1d and by using mex.

    % Interpolating all the channels at once is not much faster than doing it 1 by 1.
    F = griddedInterpolant(arrayfun(@(x)(1:x)', LFSize(1:4), 'UniformOutput', false), LF, 'linear', 'none');
    LFOut(:, :, :, USlice, :) = reshape(F(InterpIdx(2, :)', InterpIdx(1, :)', InterpIdx(4, :)', InterpIdx(3, :)'), [DestSize, NChans]);
    
    %for (ColChan = 1:NChans)

    %F = griddedInterpolant(arrayfun(@(x)(1:x)', LFSize(1:4), 'UniformOutput', false), LF, 'linear', 'none');
    %InterpSlice = F(InterpIdx(2, :)', InterpIdx(1, :)', InterpIdx(4, :)', InterpIdx(3, :)');

    %InterpSlice = interpn(squeeze(LF(:, :, :, :, ColChan)), InterpIdx(2, :), InterpIdx(1, :), InterpIdx(4, :), InterpIdx(3, :), 'linear');
    %InterpSlice = zeros(1, length(InterpIdx));
    %InterpIdx = min(LFSize(1:4)', max(1, round(InterpIdx)));
    %InterpSlice = LF(sub2ind(LFSize, InterpIdx(2, :), InterpIdx(1, :), InterpIdx(4, :), InterpIdx(3, :), repmat(ColChan, 1, length(InterpIdx))));
    %for col = 1:length(InterpIdx)
    %    InterpSlice(col) = LF(InterpIdx(2, col), InterpIdx(1, col), InterpIdx(4, col), InterpIdx(3, col), ColChan);
    %end

    %InterpSlice = reshape(InterpSlice, LFSize);
    %LFOut(:, :, :, USlice, :) = InterpSlice;
    %end

end
