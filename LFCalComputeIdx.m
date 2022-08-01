function [InterpIdx, RectOptions] = LFCalComputeIdx(LFSize, IdxSize, CalInfo, RectOptions)
    assert(length(IdxSize) == 4, 'Incorrect dimensions for IdxSize. Should be a 4 length row vector.')

    RectOptions = LFDefaultField('RectOptions', 'RectCamIntrinsicsH', LFDefaultIntrinsics(LFSize, CalInfo));
    RectOptions = LFDefaultField('RectOptions', 'Precision', 'single');

    t_in = cast(IdxSize{1}, 'uint16');
    s_in = cast(IdxSize{2}, 'uint16');
    v_in = cast(IdxSize{3}, 'uint16');
    u_in = cast(IdxSize{4}, 'uint16');
    [tt, ss, vv, uu] = ndgrid(t_in, s_in, v_in, u_in);

    InterpIdx = [ss(:)'; tt(:)'; uu(:)'; vv(:)'; ones(size(ss(:)'))];

    [InterpIdx, RectOptions] = LFMapRectifiedToMeasured(InterpIdx, CalInfo, RectOptions);
end
