function LFMetadata = LFLytroExtractMetadata(InputFname, DecodeOptions)
    %---Read the LFP or raw file + external metadata---

    % Compatibility: for loading extracted raw / json files
    DecodeOptions = LFDefaultField('DecodeOptions', 'MetadataFnamePattern', '_metadata.json');
    DecodeOptions = LFDefaultField('DecodeOptions', 'SerialdataFnamePattern', '_private_metadata.json');

    FileExtension = InputFname(end - 2:end);

    switch (lower(FileExtension))
        case 'raw' %---Load raw light field's metadata---
            FNameBase = InputFname(1:end - 4);
            MetadataFname = [FNameBase, DecodeOptions.MetadataFnamePattern];
            SerialdataFname = [FNameBase, DecodeOptions.SerialdataFnamePattern];
            LFMetadata = LFReadMetadata(MetadataFname);
            LFMetadata.SerialData = LFReadMetadata(SerialdataFname);

        otherwise %---Load Lytro LFP's metadata---
            LFP = LFReadLFP(InputFname, true);
            LFMetadata = LFP.Metadata;
            LFMetadata.SerialData = LFP.Serials;
    end

end
