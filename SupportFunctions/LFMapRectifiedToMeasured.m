% LFMapRectifiedToMeasured - Applies a calibrated camera model to map desired samples to measured samples
%
% Usage:
%
%     InterpIdx = LFMapRectifiedToMeasured( InterpIdx, CalInfo, RectOptions )
%
% Helper function used by LFCalRectifyLF.  Based on a calibrated camera model, including distortion parameters and a
% desired intrinsic matrix, the indices of a set of desired sample is mapped to the indices of corresponding measured
% samples.
%
% Inputs :
%
%     InterpIdx : Set of desired indices, in homogeneous coordinates
%     CalInfo : Calibration info as returned by LFFindCalInfo
%     RectOptions : struct controlling the rectification process, see LFCalRectify
%
% Outputs:
%
%      InterpIdx : continuous-domain indices for interpolating from the measured light field
%
% User guide: <a href="matlab:which LFToolbox.pdf; open('LFToolbox.pdf')">LFToolbox.pdf</a>
% See also: LFCalRectifyLF

% Copyright (c) 2013-2020 Donald G. Dansereau

function [InterpIdx, RectOptions] = LFMapRectifiedToMeasured(InterpIdx, CalInfo, RectOptions)
    RectOptions = LFDefaultField('RectOptions', 'NInverse_Distortion_Terms', 3);
    RectOptions = LFDefaultField('RectOptions', 'NInverse_Distortion_Iters', 2);
    RectOptions = LFDefaultField('RectOptions', 'Inverse_Distortion_Method', 'PowerSeries');
    RectOptions = LFDefaultField('RectOptions', 'Inverse_Terms', {});

    assert(RectOptions.NInverse_Distortion_Terms >= 1 && floor(RectOptions.NInverse_Distortion_Terms) ...
        == RectOptions.NInverse_Distortion_Terms, 'Number of Distortion terms must be integer and greater than 0.')

    %---Cast the to the required precision---
    InterpIdx = cast(InterpIdx, RectOptions.Precision);

    %---Convert the index of the desired ray to a ray representation using ideal intrinsics---
    InterpIdx = RectOptions.RectCamIntrinsicsH * InterpIdx;

    NumCalTerms = length(CalInfo.EstCamDistortionV);

    %---Apply inverse lens distortion to yield the undistorted ray---
    if (NumCalTerms == 5)
        k1 = CalInfo.EstCamDistortionV(1); % r^2
        k2 = CalInfo.EstCamDistortionV(2); % r^4
        k3 = CalInfo.EstCamDistortionV(3); % r^6
        b1 = CalInfo.EstCamDistortionV(4); % decentering of lens distortion
        b2 = CalInfo.EstCamDistortionV(5); % decentering of lens distortion
    else
        warning('Wrong number of distortion parameters passed')
    end

    InterpIdx(3:4, :) = InterpIdx(3:4, :) - [b1; b2];
    % InterpIdx(3:4, :) = bsxfun(@minus, InterpIdx(3:4, :), [b1; b2]); % decentering of lens distortion

    switch RectOptions.Inverse_Distortion_Method
        case 'Iterative'
            %---Iteratively estimate the undistorted direction----
            DesiredDirection = InterpIdx(3:4, :);

            for (InverseIters = 1:RectOptions.NInverse_Distortion_Iters)
                R2 = sum(InterpIdx(3:4, :).^2); % compute radius^2 for the current estimate
                % update estimate based on inverse of distortion model
                Acc = ones(size(R2));
                Prev = R2;

                for k = 1:3
                    Acc = Acc + CalInfo.EstCamDistortionV(k) .* Prev;

                    if (k == 3)
                        % Avoid unecessary multiplication
                        break
                    end

                    Prev = Prev .* R2;
                end

                % (1 + k1 .* R2 + k2 .* R2.^2 + k3 .* R2.^3);
                InterpIdx(3:4, :) = DesiredDirection ./ Acc;
            end

            clear R2 DesiredDirection
        case 'Roots'
            Ru = sqrt(sum(InterpIdx(3:4, :).^2));
            InterpIdxU = InterpIdx(3, :);
            InterpIdxV = InterpIdx(4, :);

            MultipleRoots = false(1, size(InterpIdx, 2));

            parfor col = 1:size(InterpIdx, 2)
                Rd = roots([k3, 0, k2, 0, k1, 0, 1, -Ru(col)]);
                mask = imag(Rd) == 0;
                RealRd = real(Rd(mask));
                scaler = RealRd(end) ./ Ru(col);

                MultipleRoots(col) = length(RealRd) > 1;

                InterpIdxU(col) = InterpIdxU(col) .* scaler;
                InterpIdxV(col) = InterpIdxV(col) .* scaler;
            end
            NumMultipleRoots = sum(MultipleRoots, 'all');
            if NumMultipleRoots >= 1
                warning('%d multiple roots were found.', NumMultipleRoots);
            end

            InterpIdx(3:4, :) = [InterpIdxU; InterpIdxV];
        case 'PowerSeries'
            % Taylor series is a power series...
            R2 = sum(InterpIdx(3:4, :).^2);

            BTerms = B2([k1, k2, k3], RectOptions.NInverse_Distortion_Terms);
            RectOptions.Inverse_Terms{end + 1} = BTerms;

            Q = ones(1, size(R2, 2));
            Prev = R2;

            % Can replaced with PowerseriesQ
            for i = 2:length(BTerms)
                Q = Q + BTerms(i) .* Prev;

                if (i == length(BTerms))
                    break
                end

                Prev = Prev .* R2;
            end

            InterpIdx(3:4, :) = InterpIdx(3:4, :) .* Q;

            %(B2([k1, k2, k3], NInverse_Distortion_Terms) * [ones(size(R2)); R2; R2.^2]
            %InterpIdx(3:4, :) = InterpIdx(3:4, :) .* (B2([k1, k2, k3], RectOptions.NInverse_Distortion_Terms - 1) * acc);

        case 'NewtRaph'
            Ru = sqrt(sum(InterpIdx(3:4, :).^2));
            Rd = Ru;

            for InverseIters = 1:RectOptions.NInverse_Distortion_Iters
                %F = Rd + k1 .* Rd.^3 + k2 .* Rd.^5 + k3 .* Rd.^7 - Ru;
                %Fd = 1 + 3 .* k1 .* Rd.^2 + 5 .* k2 .* Rd.^4 + 7 .* k3 .* Rd.^6;
                % Rd = Rd - (Rd + k1 .* Rd.^3 + k2 .* Rd.^5 + k3 .* Rd.^7 - Ru) ./ (1 + 3 .* k1 .* Rd.^2 + 5 .* k2 .* Rd.^4 + 7 .* k3 .* Rd.^6);
                Rd = Rd - Froot(Rd, Ru) ./ Denominiator(Rd);
            end

            % Froot = Rd + k1 .* Rd.^3 + k2 .* Rd.^5 + k3 .* Rd.^7 - Ru;

            InterpIdx(3:4, :) = InterpIdx(3:4, :) .* (Rd ./ Ru);
            % Pretty sure the below is equivalant to the above
            % InterpIdx(3:4, :) = InterpIdx(3:4, :) ./ (1 + k1 .* R2 .^ 2 + k2 .* R2.^4 + k3 .* R2.^6);
        case 'ModifiedSecant'
            Ru = sqrt(double(sum(InterpIdx(3:4, :).^2)));
            Delta = 0.001;

            Rd_1 = Ru;

            % Froot = Rd + k1 .* Rd.^3 + k2 .* Rd.^5 + k3 .* Rd.^7 - Ru;
            % Ru cancels out in denominator
            for InverseIters = 1:RectOptions.NInverse_Distortion_Iters
                Rd_2 = Rd_1 - (Froot(Rd_1, Ru) .* Delta .* Rd_1) ./ (Froot(Rd_1 + Delta .* Rd_1, 0) - Froot(Rd_1, 0));
                Rd_1 = Rd_2;
            end

            InterpIdx(3:4, :) = InterpIdx(3:4, :) .* (Rd_1 ./ Ru);
        case 'RegularFittedPowerSeries'
            Ru_requested = double(sqrt(sum(InterpIdx(3:4, :).^2)));
            Rd_fit = linspace(0, max(Ru_requested) * 2, 10000);

            % Ru_fit = Rd_fit .* (1 + k1 .* Rd_fit.^2 + k2 .* Rd_fit.^4 + k3 .* Rd_fit.^6);
            Ru_fit = Rd_fit .* PowerSeriesQ(CalInfo.EstCamDistortionV(1:3), Rd_fit.^2);

            Rd_predictor = @(x, Ru) Ru .* PowerSeriesQ(x, Ru);
            %options.Algorithm = 'levenberg-marquardt';
            options.UseParallel = true;
            options.Display = 'iter-detailed';

            lb = [];
            ub = [];

            x = lsqcurvefit(Rd_predictor, zeros(1, RectOptions.NInverse_Distortion_Terms), Ru_fit, Rd_fit, lb, ub, options);
            disp(x);
            RectOptions.Inverse_Terms{end + 1} = x;

            tStart = cputime;
            InterpIdx(3:4, :) = InterpIdx(3:4, :) .* PowerSeriesQ(x, Ru_requested);
            fprintf('\nRunning through inverse model took %.2f CPU seconds.\n', cputime - tStart);

        case 'FittedPowerSeries'

            % assert(NumCalTerms == 4 || NumCalTerms == 6, ...
            % 'Small distortion approximation only applicable when two k terms were fitted.')

            %if (NumCalTerms == 5)
            % Perform fitting to generate 'a' terms
            Ru2_requested = double(sum(InterpIdx(3:4, :).^2));

            % Run through forward model a sample of Rd values. Hence the domain of the
            % inverse will not necessarily be the set of requested Ru values.

            % Not sure if 2 * length is necessary, avoid aliasing?
            Rd_fit = linspace(0, max(Ru2_requested) * 2, 10000);
            % Ru_fit = Rd_fit .* (1 + k1 .* Rd_fit.^2 + k2 .* Rd_fit.^4 + k3 .* Rd_fit.^6);
            Ru_fit = Rd_fit .* PowerSeriesQ(CalInfo.EstCamDistortionV(1:3), Rd_fit.^2);

            % x = [a1, a2]
            %Rd_predictor = @(x, Ru) Ru + (-Ru .* (k1 .* Ru.^2 + k2.^2 .* Ru.^4 + k1 .* Ru.^4 + k2.^2 .* Ru.^8 + ...
            %2 .* k1 .* k2 .* Ru.^6) ./ (1 + 4 .* x(1) .* Ru.^2 + 6 .* x(2) .* Ru.^4));

            Rd_predictor = @(x, Ru) Ru .* PowerSeriesQ(x, Ru.^2);
            %options.Algorithm = 'levenberg-marquardt';
            options.UseParallel = true;
            options.Display = 'iter-detailed';

            lb = [];
            ub = [];

            % May want to fix the range
            % Can use nlinfit instead as Y output is a vector
            % Effect of k3 term is accounted for by the a1 and a2 terms.
            % x = lsqnonlin(RuResidualF, [0, 0], lb, ub, options);
            x = lsqcurvefit(Rd_predictor, zeros(1, RectOptions.NInverse_Distortion_Terms), Ru_fit, Rd_fit, lb, ub, options);
            disp(x);
            RectOptions.Inverse_Terms{end + 1} = x;

            tStart = cputime;
            % Therefore run through inverse model
            % Rd = Ru_requested + (-Ru_requested .* (k1 .* Ru_requested.^2 + k2.^2 .* Ru_requested.^4 + k1 .* Ru_requested.^4 + k2.^2 .* Ru_requested.^8 ...
            % + 2 .* k1 .* k2 .* Ru_requested.^6) ./ (1 + 4 .* x(1) .* Ru_requested.^2 + 6 .* x(2) .* Ru_requested.^4));

            InterpIdx(3:4, :) = InterpIdx(3:4, :) .* PowerSeriesQ(x, Ru2_requested);
            fprintf('\nRunning through inverse model took %.2f CPU seconds.\n', cputime - tStart);
            %             elseif NumCalTerms == 4 || true
            %                 % Use inverse approximation that only has two 'k' terms
            %                 warning('Ignoring k3 in inverse rectification calculation.')
            %                 Ru = sqrt(sum(InterpIdx(3:4, :).^2));
            %                 Delta_Rd = -Ru .* (k1 .* Ru.^2 + k2.^2 .* Ru.^4 + k1 .* Ru.^4 + k2.^2 .* Ru.^8 + 2 .* k1 .* k2 .* Ru.^6) ...
            %                     ./ (1 + 4 .* k1 .* Ru.^2 + 6 .* k2 .* Ru.^4);
            %                 InterpIdx(3:4, :) = InterpIdx(3:4, :) + Delta_Rd;
            %             else
            %                 % Directly use 'a' terms to calculate inverse distortion
            %                 Ru = sqrt(sum(InterpIdx(3:4, :).^2));
            %                 Delta_Rd = -Ru .* (k1 .* Ru.^2 + k2.^2 .* Ru.^4 + k1 .* Ru.^4 + k2.^2 .* Ru.^8 + 2 .* k1 .* k2 .* Ru.^6) ...
            %                     ./ (1 + 4 .* a1 .* Ru.^2 + 6 .* a2 .* Ru.^4);
            %                 InterpIdx(3:4, :) = InterpIdx(3:4, :) + Delta_Rd;
            %             end
        otherwise
            error('Unrecognised inverse distortion method.')
    end

    InterpIdx(3:4, :) = InterpIdx(3:4, :) + [b1; b2];
    % InterpIdx(3:4, :) = bsxfun(@plus, InterpIdx(3:4, :), [b1; b2]); % decentering of lens distortion

    %---Convert the undistorted ray to the corresponding index using the calibrated intrinsics---
    % todo[optimization]: The variable InterpIdx could be precomputed and saved with the calibration
    InterpIdx = CalInfo.EstCamIntrinsicsH^ - 1 * InterpIdx;

    %---Interpolate the required values---
    InterpIdx = InterpIdx(1:4, :); % drop homogeneous coordinates

    function Froot = Froot(Rd, Ru)
        % Froot = Rd + k1 .* Rd.^3 + k2 .* Rd.^5 + k3 .* Rd.^7 - Ru;

        Froot = Rd .* PowerSeriesQ(CalInfo.EstCamDistortionV(1:3), Rd.^2);

        if ~(numel(Ru) == 1 && Ru == 0)
            Froot = Froot - Ru;
        end

        %{
        Rd2 = Rd.^2;
        Prev = Rd;

        for k = 1:3
            Prev = Prev .* Rd2;
            Froot = Froot + CalInfo.EstCamDistortionV(k) .* Prev;
        end
        %}
    end

    function Denominiator = Denominiator(Rd)
        % Denominiator = 1 + 3 .* k1 .* Rd.^2 + 5 .* k2 .* Rd.^4 + 7 .* k3 .* Rd.^6
        Derivative = [3, 5, 7];

        Denominiator = PowerSeriesQ(CalInfo.EstCamDistortionV(1:3) .* Derivative, Rd.^2);

        %{
        Rd2 = Rd.^2;
        Prev = Rd2;
        Denominiator = ones(size(Rd));

        for k = 1:3
            Denominiator = Denominiator + Derivative(k) .* CalInfo.EstCamDistortionV(k) .* Rd2;
            if k == 3
                % Avoid unnecessary multiplication
                break
            end
            Prev = Prev .* Rd2;
        end
        %}
    end

end

function Q = PowerSeriesQ(x, Ru)
    % Q = 1 + k1 .* Ru. + k2 .* Ru.^2 ...
    Q = ones(size(Ru));
    Prev = Ru;

    for i = 1:size(x, 2)
        Q = Q + x(i) .* Prev;

        if (i == size(x, 2))
            % Avoid unnecessary multiplication
            break
        end

        Prev = Prev .* Ru;
    end

end

function p = B1(a, N)
    % p = zeros(N, 2 * N - 1);

    function return_val = p_get(j, k)

        if j >= 2 && k <= 1 || j < 1
            return_val = 0;
        else
            return_val = p(j, k);
        end

    end

    for k = 1:2 * N - 1
        p(1, k) = 1;
    end

    for k = 2:2 * N - 1

        for j = 2:N
            %p(j, k) = p_get(j, k - 1) + a(1) * p_get(j - 1, k - 1) + a(2) * ...
            %    p_get(j - 2, k - 1) + a(3) * p_get(j - 3, k - 1) ...
            %    + a(4) * p_get(j - 4, k - 1);

            p(j, k) = p_get(j, k - 1);

            for i = 1:length(a)
                p(j, k) = p(j, k) + a(i) * p_get(j - i, k - 1);
            end

        end

    end

end

function b = B2(a, N)
    q = zeros(1, length(a));
    q(1) = 1;
    b = zeros(1, N + 1);

    b(1) = 1;

    if (N == 0)
        return;
    end

    p = B1(a, N);

    for n = 1:N
        tmp = 0;

        for k = 1:length(a)
            tmp = tmp - a(k) * q(k);
        end

        b(n + 1) = tmp;

        for j = 1:8/9*n
            
            k = n - j;
            b(n + 1) = b(n + 1) - b(k + 1) * p(n - k + 1, 2 * k + 1);

        end

        q(2:length(a)) = q(1:length(a) - 1);
        q(1) = tmp;
    end

end
