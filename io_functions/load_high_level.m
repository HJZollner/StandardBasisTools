function [BasisMatrix, ProvenanceJson] = load_high_level(HighLevelData,JsonFile)
%% [BasisMatrix] = load_high_level(HighLevelData,JsonFile)
%   This function loads the high-level basis set file 
%
%   USAGE:
%       [BasisMatrix] = load_high_level(HighLevelData,JsonFile)
%
%   INPUTS:
%       HighLevelData   = input data in matlab structure format.
%       JsonFile        = JsonFile with provanance 
%
%   OUTPUTS:
%       BasisMatrix     = basis matrix
%       ProvenanceJson  = simualtion provenance json
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
if nargin < 2
    error('You need to specifiy a high-level basis file and a simualtion provenance file!');
end

%% Load data using the npy-matlab package and standard json decoding
BasisMatrix = readNPY(HighLevelData);                                       % Load basis matrix from standardized npy file
                                                                            % Load json and store provenance


end

