%% Minimum example to generate standardized low-level basis set and LCM fitting results outputs
% 
% PREREQUISITS:  A recent verison of Osprey (https://github.com/schorschinho/osprey) from the develop branch
%% Analyze example HERCULES data 
%PLEASE CHANGE PATH ACCORDING TO YOUR GIT
pathToData = '/Users/helge/Documents/GitHub/StandardBasisTools/StandardBasisTools/example/NIfTi';


MRSCont = OspreyJob(fullfile(pathToData,'jobNIfTIMRSHERC.json'),0,'11'); % Initialize jobfile
MRSCont = OspreyLoad(MRSCont); % Load data
MRSCont = OspreyProcess(MRSCont); % Process data

% Add processing provenance
% Metabolite data
hdr_ext = MRSCont.processed.metab{1, 1}.nii_mrs.hdr_ext;
hdr_ext.ProcessingApplied(1).Time = datestr(now,30);
hdr_ext.ProcessingApplied(1).Program = 'Osprey';
hdr_ext.ProcessingApplied(1).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(1).Method = 'RF coil combination';
hdr_ext.ProcessingApplied(1).Details = 'Weighted summation (Hall et al. 2014), reference = raw_ref';

hdr_ext.ProcessingApplied(2).Time = datestr(now,30);
hdr_ext.ProcessingApplied(2).Program = 'Osprey';
hdr_ext.ProcessingApplied(2).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(2).Method = 'Frequency and phase correction';
hdr_ext.ProcessingApplied(2).Details = 'Robust Spectral Registration (Mikkelsen et al. 2020), dim =  DIM_DYN';

hdr_ext.ProcessingApplied(3).Time = datestr(now,30);
hdr_ext.ProcessingApplied(3).Program = 'Osprey';
hdr_ext.ProcessingApplied(3).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(3).Method = 'Signal averaging';
hdr_ext.ProcessingApplied(3).Details = 'Similarity metric weighted sum (Mikkelsen et al. 2020)';

hdr_ext.ProcessingApplied(4).Time = datestr(now,30);
hdr_ext.ProcessingApplied(4).Program = 'Osprey';
hdr_ext.ProcessingApplied(4).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(4).Method = 'Eddy current correction';
hdr_ext.ProcessingApplied(4).Details = 'ECC (Klose 1990), reference = raw_ref';

hdr_ext.ProcessingApplied(5).Time = datestr(now,30);
hdr_ext.ProcessingApplied(5).Program = 'Osprey';
hdr_ext.ProcessingApplied(5).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(5).Method = 'Aligment of subtraction sub-spectra';
hdr_ext.ProcessingApplied(5).Details = 'L1 norm optimization of HADMARD spectra';

hdr_ext.ProcessingApplied(6).Time = datestr(now,30);
hdr_ext.ProcessingApplied(6).Program = 'Osprey';
hdr_ext.ProcessingApplied(6).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(6).Method = 'Nuisance peak removal';
hdr_ext.ProcessingApplied(6).Details = 'HSVD removal of residual water (Barkhuijsen et al 1987)';

hdr_ext.ProcessingApplied(7).Time = datestr(now,30);
hdr_ext.ProcessingApplied(7).Program = 'Osprey';
hdr_ext.ProcessingApplied(7).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(7).Method = 'Addition of sub-spectra';
hdr_ext.ProcessingApplied(7).Details = 'HERCULES Hadmard combination (Oeltzschner et al. 2019)';
MRSCont.processed.metab{1, 1}.nii_mrs.hdr_ext = hdr_ext;

% Add processing provenance
% Water reference data
hdr_ext = MRSCont.processed.ref{1, 1}.nii_mrs.hdr_ext;
hdr_ext.ProcessingApplied(1).Time = datestr(now,30);
hdr_ext.ProcessingApplied(1).Program = 'Osprey';
hdr_ext.ProcessingApplied(1).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(1).Method = 'RF coil combination';
hdr_ext.ProcessingApplied(1).Details = 'Weighted summation (Hall et al. 2014), reference = raw_ref';

hdr_ext.ProcessingApplied(2).Time = datestr(now,30);
hdr_ext.ProcessingApplied(2).Program = 'Osprey';
hdr_ext.ProcessingApplied(2).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(2).Method = 'Frequency and phase correction';
hdr_ext.ProcessingApplied(2).Details = 'Spectral Registration (Near et al. 2015), dim =  DIM_DYN';

hdr_ext.ProcessingApplied(3).Time = datestr(now,30);
hdr_ext.ProcessingApplied(3).Program = 'Osprey';
hdr_ext.ProcessingApplied(3).Version = MRSCont.ver.Osp;
hdr_ext.ProcessingApplied(3).Method = 'Signal averaging';
hdr_ext.ProcessingApplied(3).Details = 'Simple signal averaging';
MRSCont.processed.ref{1, 1}.nii_mrs.hdr_ext = hdr_ext;
[MRSCont] = osp_saveNII(MRSCont);


MRSCont = OspreyFit(MRSCont); %Fit data
%% Prepare NIfTI-MRS fitting output
% This is a preliminary implementation of the proposed guidelines
% Gather metabolite data fits from GABA-edited and GSH-edited difference and sum spectra
which_spec = 'metab';
for sub = 1 : 3
    VoxelIndex = [1 sub 1];
    dataToPlot=MRSCont.processed.metab{1};
    fitRangePPM = MRSCont.opts.fit.range;
    basisSet    = MRSCont.fit.resBasisSet.([which_spec]).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]){VoxelIndex(3),1,VoxelIndex(2)};
    dataToPlot   = op_takesubspec(MRSCont.processed.(which_spec){1},find(strcmp(MRSCont.processed.(which_spec){1}.names,basisSet.names{1})));
     fitParams   = MRSCont.fit.results.(which_spec).fitParams{VoxelIndex(3),1,VoxelIndex(2)};
    inputData.dataToFit                 = dataToPlot;
    inputData.basisSet                  = basisSet;
    inputSettings.scale                 = MRSCont.fit.scale{1};
    inputSettings.fitRangePPM           = fitRangePPM;
    inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
    inputSettings.fitStyle              = MRSCont.opts.fit.style;
    inputSettings.flags.isMEGA          = MRSCont.flags.isMEGA;
    inputSettings.flags.isHERMES        = MRSCont.flags.isHERMES;
    inputSettings.flags.isHERCULES      = MRSCont.flags.isHERCULES;
    inputSettings.flags.isPRIAM         = MRSCont.flags.isPRIAM;
    inputSettings.concatenated.Subspec  = 'diff1';
    if isfield(MRSCont.opts.fit,'GAP') && ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
        inputSettings.GAP = MRSCont.opts.fit.GAP.(dataToPlot.names{1});
    else
        inputSettings.GAP = [];
    end
    [ModelOutput{sub}] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams);
end

% For testing purposes only
% figure, plot(ModelOutput{1}.ppm,ModelOutput{1}.data) 

% Package NIfTI results from the metabolite fit
NIfTIfits = MRSCont.processed.metab{1, 1};
NIfTIfits.specs =[];
NIfTIfits.specs(:,1,1) = ModelOutput{1}.data;
NIfTIfits.specs(:,2,1) = ModelOutput{1}.completeFit;
NIfTIfits.specs(:,3,1) = ModelOutput{1}.residual;
NIfTIfits.specs(:,4,1) = ModelOutput{1}.baseline;
NIfTIfits.specs(:,5:5+size(ModelOutput{1}.indivMets,2)-1,1) = ModelOutput{1}.indivMets(:,:);
NIfTIfits.ppm = ModelOutput{1}.ppm;

NIfTIfits.specs(:,1,2) = ModelOutput{2}.data;
NIfTIfits.specs(:,2,2) = ModelOutput{2}.completeFit;
NIfTIfits.specs(:,3,2) = ModelOutput{2}.residual;
NIfTIfits.specs(:,4,2) = ModelOutput{2}.baseline;
NIfTIfits.specs(:,5:5+size(ModelOutput{2}.indivMets,2)-1,2) = ModelOutput{2}.indivMets(:,:);

NIfTIfits.specs(:,1,3) = ModelOutput{3}.data;
NIfTIfits.specs(:,2,3) = ModelOutput{3}.completeFit;
NIfTIfits.specs(:,3,3) = ModelOutput{3}.residual;
NIfTIfits.specs(:,4,3) = ModelOutput{3}.baseline;
NIfTIfits.specs(:,5:5+size(ModelOutput{3}.indivMets,2)-1,3) = ModelOutput{3}.indivMets(:,:);

NIfTIfits.fids = ifft(fftshift(NIfTIfits.specs,1),[],1);
NIfTIfits.specs = fftshift(fft(NIfTIfits.fids,[],1),1);
NIfTIfits.sz = size(NIfTIfits.specs);
f = (NIfTIfits.ppm) * NIfTIfits.Bo*42.577;
NIfTIfits.spectralwidth= f(end)-f(1);
dt = 1/NIfTIfits.spectralwidth;
NIfTIfits.dwelltime = 1/NIfTIfits.spectralwidth;
NIfTIfits.t   = [0 : dt : (NIfTIfits.sz(1)-1)*dt];
NIfTIfits.nii_mrs.hdr_ext.SpectralDwellTime =NIfTIfits.dwelltime;
NIfTIfits.nii_mrs.hdr.pixdim(5) = NIfTIfits.dwelltime;
NIfTIfits.dims.subSpecs=3;
NIfTIfits.dims.averages=2;
NIfTIfits.nii_mrs.hdr_ext =rmfield(NIfTIfits.nii_mrs.hdr_ext,'dim_5');
NIfTIfits.nii_mrs.hdr_ext =rmfield(NIfTIfits.nii_mrs.hdr_ext,'dim_6');
NIfTIfits.nii_mrs.hdr_ext =rmfield(NIfTIfits.nii_mrs.hdr_ext,'dim_7');

hdr_ext = NIfTIfits.nii_mrs.hdr_ext;
hdr_ext.FittingPerformed(1).Time = datestr(now,30);
hdr_ext.FittingPerformed(1).Program = 'Osprey';
hdr_ext.FittingPerformed(1).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(1).Method = 'Zero-filling';
hdr_ext.FittingPerformed(1).Details = 'Zerofill twice';

basisNames = MRSCont.fit.resBasisSet.metab.np_sw_4096_4000{1}.name;
hdr_ext.FittingPerformed(2).Time = datestr(now,30);
hdr_ext.FittingPerformed(2).Program = 'Osprey';
hdr_ext.FittingPerformed(2).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(2).Method = 'Model components DIM_EDIT 1';
hdr_ext.FittingPerformed(2).Details = cell2mat(strcat(basisNames, [mat2cell(repmat(',',1,length(basisNames)-1), 1, ones(1,length(basisNames)-1)),{''}]));

basisNames = MRSCont.fit.resBasisSet.metab.np_sw_4096_4000{2}.name;
hdr_ext.FittingPerformed(3).Time = datestr(now,30);
hdr_ext.FittingPerformed(3).Program = 'Osprey';
hdr_ext.FittingPerformed(3).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(3).Method = 'Model components DIM_EDIT 2';
hdr_ext.FittingPerformed(3).Details = cell2mat(strcat(basisNames, [mat2cell(repmat(',',1,length(basisNames)-1), 1, ones(1,length(basisNames)-1)),{''}]));

basisNames = MRSCont.fit.resBasisSet.metab.np_sw_4096_4000{3}.name;
hdr_ext.FittingPerformed(4).Time = datestr(now,30);
hdr_ext.FittingPerformed(4).Program = 'Osprey';
hdr_ext.FittingPerformed(4).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(4).Method = 'Model components DIM_EDIT 3';
hdr_ext.FittingPerformed(4).Details = cell2mat(strcat(basisNames, [mat2cell(repmat(',',1,length(basisNames)-1), 1, ones(1,length(basisNames)-1)),{''}]));

hdr_ext.FittingPerformed(5).Time = datestr(now,30);
hdr_ext.FittingPerformed(5).Program = 'Osprey';
hdr_ext.FittingPerformed(5).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(5).Method = 'Physics Model';
hdr_ext.FittingPerformed(5).Details = 'Default Osprey, Spectral Domain Optimization with unregularized spline baseline';

hdr_ext.FittingPerformed(6).Time = datestr(now,30);
hdr_ext.FittingPerformed(6).Program = 'Osprey';
hdr_ext.FittingPerformed(6).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(6).Method = 'LCM settings';
hdr_ext.FittingPerformed(6).Details = ['Default Osprey LCM settings bspl = ' MRSCont.opts.fit.bLineKnotSpace ', Style = ' MRSCont.opts.fit.style];

hdr_ext.FittingPerformed(7).Time = datestr(now,30);
hdr_ext.FittingPerformed(7).Program = 'Osprey';
hdr_ext.FittingPerformed(7).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(7).Method = 'Basis set file';
[~,basis,~] = fileparts(MRSCont.opts.fit.basisSetFile);
hdr_ext.FittingPerformed(5).Details = ['Basis Set File ' basis];
NIfTIfits.nii_mrs.hdr_ext = hdr_ext;

nii = io_writeniimrs(NIfTIfits, fullfile(pathToData,'derivatives/NIfTIMRS/sub_01_ses-01_mrs_sub_01_ses-01_acq-hercules_type-act_svs-fit.nii.gz'),{'DIM_EDIT'});
% For testing purposes only
% figure, plot(NIfTIfits.ppm,squeeze(NIfTIfits.specs(:,1,1))), hold on, plot(NIfTIfits.ppm,squeeze(NIfTIfits.specs(:,1,2))), plot(NIfTIfits.ppm,squeeze(NIfTIfits.specs(:,1,3))); 
% nii_fits_load = io_loadspec_niimrs(fullfile(pathToData,'derivatives/NIfTIMRS/sub_01_ses-01_mrs_sub_01_ses-01_acq-hercules_type-act_svs-fit.nii.gz'));
% figure, plot(nii_fits_load.ppm-2.48,squeeze(nii_fits_load.specs(:,1,1))), hold on, plot(nii_fits_load.ppm-2.48,squeeze(nii_fits_load.specs(:,1,2))),plot(nii_fits_load.ppm-2.48,squeeze(nii_fits_load.specs(:,1,3)));



%% Gather water data fits
which_spec = 'ref';
VoxelIndex = [1 1 1];
dataToPlot=MRSCont.processed.metab{1};
fitRangePPM = MRSCont.opts.fit.range;
basisSet    = MRSCont.fit.resBasisSet.([which_spec]).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]){VoxelIndex(3),1,VoxelIndex(2)};
dataToPlot   = op_takesubspec(MRSCont.processed.(which_spec){1},find(strcmp(MRSCont.processed.(which_spec){1}.names,basisSet.names{1})));
 fitParams   = MRSCont.fit.results.(which_spec).fitParams{VoxelIndex(3),1,VoxelIndex(2)};
inputData.dataToFit                 = dataToPlot;
inputData.basisSet                  = basisSet;
inputSettings.scale                 = MRSCont.fit.scale{1};
inputSettings.fitRangePPM           = fitRangePPM;
inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
inputSettings.fitStyle              = MRSCont.opts.fit.style;
inputSettings.flags.isMEGA          = MRSCont.flags.isMEGA;
inputSettings.flags.isHERMES        = MRSCont.flags.isHERMES;
inputSettings.flags.isHERCULES      = MRSCont.flags.isHERCULES;
inputSettings.flags.isPRIAM         = MRSCont.flags.isPRIAM;
inputSettings.concatenated.Subspec  = 'diff1';
if isfield(MRSCont.opts.fit,'GAP') && ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
    inputSettings.GAP = MRSCont.opts.fit.GAP.(dataToPlot.names{1});
else
    inputSettings.GAP = [];
end
[ModelOutput{1}] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
% figure, plot(ModelOutput{1}.ppm,ModelOutput{1}.data) 

% Package NIfTI results from the metabolite fit
NIfTIfits = MRSCont.processed.ref{1, 1};
NIfTIfits.specs =[];
NIfTIfits.specs(:,1) = ModelOutput{1}.data;
NIfTIfits.specs(:,2) = ModelOutput{1}.completeFit;
NIfTIfits.specs(:,3) = ModelOutput{1}.residual;
NIfTIfits.ppm = ModelOutput{1}.ppm;


NIfTIfits.fids = ifft(fftshift(NIfTIfits.specs,1),[],1);
NIfTIfits.specs = fftshift(fft(NIfTIfits.fids,[],1),1);
NIfTIfits.sz = size(NIfTIfits.specs);
f = (NIfTIfits.ppm) * NIfTIfits.Bo*42.577;
NIfTIfits.spectralwidth= f(end)-f(1);
dt = 1/NIfTIfits.spectralwidth;
NIfTIfits.dwelltime = 1/NIfTIfits.spectralwidth;
NIfTIfits.t   = [0 : dt : (NIfTIfits.sz(1)-1)*dt];
NIfTIfits.nii_mrs.hdr_ext.SpectralDwellTime =NIfTIfits.dwelltime;
NIfTIfits.nii_mrs.hdr.pixdim(5) = NIfTIfits.dwelltime;
NIfTIfits.dims.averages=2;
NIfTIfits.nii_mrs.hdr_ext =rmfield(NIfTIfits.nii_mrs.hdr_ext,'dim_5');
NIfTIfits.nii_mrs.hdr_ext =rmfield(NIfTIfits.nii_mrs.hdr_ext,'dim_6');
NIfTIfits.nii_mrs.hdr_ext =rmfield(NIfTIfits.nii_mrs.hdr_ext,'dim_7');

hdr_ext = NIfTIfits.nii_mrs.hdr_ext;
hdr_ext.FittingPerformed(1).Time = datestr(now,30);
hdr_ext.FittingPerformed(1).Program = 'Osprey';
hdr_ext.FittingPerformed(1).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(1).Method = 'Zero-filling';
hdr_ext.FittingPerformed(1).Details = 'Zerofill twice';

basisNames = MRSCont.fit.resBasisSet.ref.np_sw_4096_4000{1}.name;
hdr_ext.FittingPerformed(2).Time = datestr(now,30);
hdr_ext.FittingPerformed(2).Program = 'Osprey';
hdr_ext.FittingPerformed(2).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(2).Method = 'Model components';
hdr_ext.FittingPerformed(2).Details = cell2mat(strcat(basisNames, [mat2cell(repmat(',',1,length(basisNames)-1), 1, ones(1,length(basisNames)-1)),{''}]));



hdr_ext.FittingPerformed(5).Time = datestr(now,30);
hdr_ext.FittingPerformed(5).Program = 'Osprey';
hdr_ext.FittingPerformed(5).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(5).Method = 'Physics Model';
hdr_ext.FittingPerformed(5).Details = 'Default Osprey, Spectral Domain water fit';

hdr_ext.FittingPerformed(6).Time = datestr(now,30);
hdr_ext.FittingPerformed(6).Program = 'Osprey';
hdr_ext.FittingPerformed(6).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(6).Method = 'LCM settings';
hdr_ext.FittingPerformed(6).Details = ['Default Osprey LCM water fit'];

hdr_ext.FittingPerformed(7).Time = datestr(now,30);
hdr_ext.FittingPerformed(7).Program = 'Osprey';
hdr_ext.FittingPerformed(7).Version = MRSCont.ver.Osp;
hdr_ext.FittingPerformed(7).Method = 'Basis set file';
[~,basis,~] = fileparts(MRSCont.opts.fit.basisSetFile);
hdr_ext.FittingPerformed(5).Details = ['Basis Set File ' basis];
NIfTIfits.nii_mrs.hdr_ext = hdr_ext;

nii = io_writeniimrs(NIfTIfits, fullfile(pathToData,'/derivatives/NIfTIMRS/sub_01_ses-01_mrs_sub_01_ses-01_acq-hercules_type-ref_svs_REF-fit.nii.gz'));
% For testing purposes only
% figure, plot(NIfTIfits.ppm,squeeze(NIfTIfits.specs(:,1,1))), hold on, plot(NIfTIfits.ppm,squeeze(NIfTIfits.specs(:,1,2))), plot(NIfTIfits.ppm,squeeze(NIfTIfits.specs(:,1,3))); 
% nii_load = io_loadspec_niimrs(fullfile(pathToData,'/derivatives/NIfTIMRS/sub_01_ses-01_mrs_sub_01_ses-01_acq-hercules_type-ref_svs_REF-fit.nii.gz'));
% figure, plot(nii_load.ppm-2.48,squeeze(nii_load.specs(:,1,1))), hold on, plot(nii_load.ppm-2.48,squeeze(nii_load.specs(:,1,2))),plot(nii_load.ppm-2.48,squeeze(nii_load.specs(:,1,3)));


%% Write low-level basis spectra
% We are using the reampled basis spectra from Osprey here


NIfTIbasis = MRSCont.processed.metab{1, 1};
NIfTIbasis.specs =[];
basiset = MRSCont.fit.resBasisSet.metab.np_sw_4096_4000{1};
NIfTIbasis.specs(:,1:25,1) = basiset.specs;
basiset = MRSCont.fit.resBasisSet.metab.np_sw_4096_4000{2};
NIfTIbasis.specs(:,1:26,2) = basiset.specs;
basiset = MRSCont.fit.resBasisSet.metab.np_sw_4096_4000{3};
NIfTIbasis.specs(:,1:25,3) = basiset.specs;

NIfTIbasis.ppm = ModelOutput{1}.ppm;


NIfTIbasis.fids = ifft(fftshift(NIfTIbasis.specs,1),[],1);
NIfTIbasis.sz = size(NIfTIbasis.specs);
NIfTIbasis.spectralwidth= basiset.spectralwidth;
NIfTIbasis.dwelltime = basiset.dwelltime;
NIfTIbasis.t   = basiset.t;
NIfTIbasis.nii_mrs.hdr_ext.SpectralDwellTime =NIfTIbasis.dwelltime;
NIfTIbasis.nii_mrs.hdr.pixdim(5) = NIfTIbasis.dwelltime;
NIfTIbasis.dims.averages=2;
NIfTIbasis.dims.subSpecs=3;
NIfTIbasis.nii_mrs.hdr_ext =rmfield(NIfTIbasis.nii_mrs.hdr_ext,'dim_5');
NIfTIbasis.nii_mrs.hdr_ext =rmfield(NIfTIbasis.nii_mrs.hdr_ext,'dim_6');
NIfTIbasis.nii_mrs.hdr_ext =rmfield(NIfTIbasis.nii_mrs.hdr_ext,'dim_7');
NIfTIbasis.nii_mrs.hdr_ext =rmfield(NIfTIbasis.nii_mrs.hdr_ext,'ProcessingApplied');

hdr_ext = NIfTIbasis.nii_mrs.hdr_ext;
hdr_ext.BasisSetSimulation.Steps.Time = datestr(now,30);
hdr_ext.BasisSetSimulation.Steps.Program = 'MRS Cloud';
hdr_ext.BasisSetSimulation.Steps.Version = '1.0';
hdr_ext.BasisSetSimulation.Steps.Method = 'Fully localized matrix density simulation';
hdr_ext.BasisSetSimulation.Steps.Details = 'Cloud implementation of FID-A (Simpson et al. 2015)';

basisNames = MRSCont.fit.resBasisSet.metab.np_sw_4096_4000{1}.name;
hdr_ext.BasisSetSimulation.Metabs{1} = cell2mat(strcat(basisNames, [mat2cell(repmat(',',1,length(basisNames)-1), 1, ones(1,length(basisNames)-1)),{''}]));
basisNames = MRSCont.fit.resBasisSet.metab.np_sw_4096_4000{2}.name;
hdr_ext.BasisSetSimulation.Metabs{2} = cell2mat(strcat(basisNames, [mat2cell(repmat(',',1,length(basisNames)-1), 1, ones(1,length(basisNames)-1)),{''}]));
basisNames = MRSCont.fit.resBasisSet.metab.np_sw_4096_4000{3}.name;
hdr_ext.BasisSetSimulation.Metabs{3} = cell2mat(strcat(basisNames, [mat2cell(repmat(',',1,length(basisNames)-1), 1, ones(1,length(basisNames)-1)),{''}]));

hdr_ext.BasisSetSimulation.RefPulse.PulseShape = 'univ_eddenrefo.pta';
hdr_ext.BasisSetSimulation.RefPulse.PulseDuration = 0.007;
hdr_ext.BasisSetSimulation.RefPulse.Nucleus = '1H';

hdr_ext.BasisSetSimulation.EditPulse{1}.PulseShape = 'sl_univ_pulse.pta';
hdr_ext.BasisSetSimulation.EditPulse{1}.PulseOffset = 4.58;
hdr_ext.BasisSetSimulation.EditPulse{1}.PulseDuration = 0.02;
hdr_ext.BasisSetSimulation.EditPulse{1}.Nucleus = '1H';

hdr_ext.BasisSetSimulation.EditPulse{2}.PulseShape = 'sl_univ_pulse.pta';
hdr_ext.BasisSetSimulation.EditPulse{2}.PulseOffset = 4.18;
hdr_ext.BasisSetSimulation.EditPulse{2}.PulseDuration = 0.02;
hdr_ext.BasisSetSimulation.EditPulse{2}.Nucleus = '1H';

hdr_ext.BasisSetSimulation.EditPulse{3}.PulseShape = 'dl_Philips_4_58_1_90.pta';
hdr_ext.BasisSetSimulation.EditPulse{3}.PulseOffset = [4.58, 1.9];
hdr_ext.BasisSetSimulation.EditPulse{3}.PulseDuration = 0.02;
hdr_ext.BasisSetSimulation.EditPulse{3}.Nucleus = '1H';

hdr_ext.BasisSetSimulation.EditPulse{4}.PulseShape = 'dl_Philips_4_18_1_90.pta';
hdr_ext.BasisSetSimulation.EditPulse{4}.PulseOffset = [4.18, 1.9];
hdr_ext.BasisSetSimulation.EditPulse{4}.PulseDuration = 0.02;
hdr_ext.BasisSetSimulation.EditPulse{4}.Nucleus = '1H';

hdr_ext.BasisSetSimulation.Timing.TE1 = 0.00655*2;
hdr_ext.BasisSetSimulation.Timing.TE2 = 0.0669;


NIfTIbasis.nii_mrs.hdr_ext = hdr_ext;

nii = io_writeniimrs(NIfTIbasis, fullfile(pathToData,'/HERCULES_basis.nii.gz'),{'DIM_EDIT'});
% For testing purposes only
% figure, plot(NIfTIfits.ppm,squeeze(NIfTIfits.specs(:,1,1))), hold on, plot(NIfTIfits.ppm,squeeze(NIfTIfits.specs(:,1,2))), plot(NIfTIfits.ppm,squeeze(NIfTIfits.specs(:,1,3))); 
% nii_load = io_loadspec_niimrs(fullfile(pathToData,'/HERCULES_basis.nii.gz'));
% figure, plot(nii_load.ppm-2.48,squeeze(nii_load.specs(:,1,1))), hold on, plot(nii_load.ppm-2.48,squeeze(nii_load.specs(:,1,2))),plot(nii_load.ppm-2.48,squeeze(nii_load.specs(:,1,3)));