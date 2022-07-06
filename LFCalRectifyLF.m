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

function [LF, RectOptions] = LFCalRectifyLF(LF, CalInfo, RectOptions)

    %---Defaults---
    RectOptions = LFDefaultField('RectOptions', 'MaxUBlkSize', 32);
    RectOptions = LFDefaultField('RectOptions', 'Precision', 'single');
    RectOptions = LFDefaultField('RectOptions', 'CachedIdx', []);
    LFSize = size(LF);

    %---Build interpolation indices---
    fprintf('Generating interpolation indices...\n');
    LF = cast(LF, RectOptions.Precision);

    fprintf('Interpolating...');

    LFOut = LF;
    
    if ~isempty(RectOptions.CachedIdx)
        % Use the cached index values

        % Interpolate the entire light field in one go
        LFOut = InterpolateColours(LF, LFOut, RectOptions.CachedIdx, LFSize(1:4), 1:LFSize(4));

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
            LFOut = InterpolateColours(LF, LFOut, LFCalComputeIdx(LFSize, [arrayfun(@(x)1:x, LFSize(1:3), 'UniformOutput', false), ... 
                UStart:UStop], CalInfo, RectOptions), DestSize, UStart:UStop);

            fprintf('.')
        end

    end

    LF = LFOut;
    clear LFOut;

    %---Clip interpolation result, which sometimes rings slightly out of range---
    LF(isnan(LF)) = 0;
    LF = max(0, min(1, LF));

    fprintf('\nDone\n');
end

function LFOut = InterpolateColours(LF, LFOut, InterpIdx, DestSize, USlice)
    LFSize = size(LF);
    NChans = LFSize(5);

    for (ColChan = 1:NChans)
        % todo[optimization]: use a weighted interpolation scheme to exploit the weight channel
        InterpSlice = interpn(squeeze(LF(:, :, :, :, ColChan)), InterpIdx(2, :), InterpIdx(1, :), InterpIdx(4, :), InterpIdx(3, :), 'linear');
        InterpSlice = reshape(InterpSlice, DestSize);
        LFOut(:, :, :, USlice, ColChan) = InterpSlice;
    end

end
