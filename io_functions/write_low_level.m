function [basis_struct] = write_low_level(BasisMatrix,ProvenanceJson,SpectralWidth,DataPoints,ToExport,LB)
%% [basis_struct] = write_low_level(BasisMatrix,ProvenanceJson)
%   This function converts the high-level basis set file and converts it
%   to tool-specific low-level basis files
%
%   USAGE:
%       [BasisMatrix] = load_high_level(HighLevelData,JsonFile)
%
%   INPUTS:
%       BasisMatrix     = basis matrix
%       ProvenanceJson  = simualtion provenance json
%       SpectralWidth   = target spectral width
%       DataPoints      = target data points
%       ToExport        = file formats to export
%       LB              = apodization applied
%
%   OUTPUTS:
%       BasisMatrix     = basis matrix
%       ProvenanceJson  = simulation provenance json
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-05-31)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%       This code is also based on numerous functions from the spec2nii
%       conversion tool by Dr. William Clarke (University of Oxford)
%       https://github.com/wtclarke/spec2nii
%       Clarke et al., Magn Reson Med 88:2358-2370 (2022)
%
%   HISTORY:
%       2023-05-31: First version of the code.
%% Error handling
if nargin < 6
    LB = 3;                                                                 % Apply 3 Hz LB by default
    if nargin < 5
        ToExport.LCModel    = 1;
        ToExport.Tarquin    = 1;
        ToExport.jMRUI      = 1;
        ToExport.INSPECTOR  = 1;
        ToExport.VIFF       = 0;
        ToExport.Osprey     = 1;
        ToExport.NIfTI      = 1;
        if nargin < 4
            DataPoints = 2048;
            if nargin < 3
                     SpectralWidth =  2000;                                             
                if nargin < 2
                    error('Run the load_high_level function first!');
                end
            end
        end    
    end
end
%% Convert data to FIDs if needed
%% Apply line broadening 
%% Export calls go here
end

