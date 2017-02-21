function [success surfTimes] = get_coh_task(param)
% [success surfTimes] = get_heights_task(param)
%
% Cluster task for get_heights. Does the actual data loading
% and surface tracking.
%
% param = struct controlling the loading, processing, surface tracking,
%   and quick look generation
%  .load = structure for which records to load
%   .records_fn = filename of records file
%   .recs = current records
%   .imgs = cell vector of images to load, each image is Nx2 array of
%     wf/adc pairs
%     NOTE: wfs/adc pairs are not indices into anything, they are the absolute
%     waveform/adc numbers. The records file will be loaded to decode
%     which index each wf/adc belongs to.
%  .debug_level = debug level (scalar integer)
%
%  .proc = structure containing information about framing
%   .frm = only used to determine the filename of the output
%
% .get_heights = structure controlling get_heights processing
%  .radar_name = name of radar string
%  .season_name = name of mission string
%  .day_seg = day-segment string
%  
%  get_heights fields used by load_mcords_data.m (see that function for details)
%  .ft_wind
%  .ft_wind_time
%  .trim_vals
%  .pulse_rfi.en
%  .pulse_rfi.inc_ave
%  .pulse_rfi.thresh_scale
%  .radar
%   .Vpp_scale = radar Vpp for full scale quanitization
%   .rxs = struct array of receiver equalization coefficients
%    .wfs = struct array for each waveform
%     .chan_equal = scalar complex double (data DIVIDED by this)
%     .td = time delay correction (applied during pulse compression)
%  
%  get_heights fields for post processing
%  .roll_correction = boolean, whether or not to apply roll phase correction
%  .lever_arm_fh = string containing function name
%  .elev_correction = boolean, whether or not to apply Elevation1 phase correction
%  .B_filter = double vector, FIR filter coefficients to apply before
%    decimating, this function loads data before and after this frame
%    (if available) to avoid transients at the beginning and end
%  .decimate_factor = positive integer, decimation rate
%  .inc_B_filter: double vector, FIR filter coefficients to apply before
%    incoherent average decimation. If not defined or empty, then
%    inc_B_filter is set to ones(1,inc_ave)/inc_ave.
%  .inc_ave = integer scalar, number of incoherent averages to apply
%    (also decimates by this number). If set to < 1, complex data are
%    returned.  Setting to 1 causes the data to be power detected (i.e.
%    become incoherent), but no averaging is done.
%
%  .surf = get_heights structure controlling surface tracking
%   .en = boolean, whether or not to apply surface tracking
%   .wf_idx = positive integer, which waveform in the wfs list to apply surface tracking on
%   .min_bin = double scalar, the minimum range time that the surface can be tracked to.
%     This is used to keep the surface tracking routine from picking up the
%     feedthrough.  It requires a minimum Elevation1 AGL.
%   .manual = boolean, whether or not to enable the manual tracking
%     interface.  Generally better to let the automated routine run, fix in
%     picker, and then update records (so surf.manual is mostly for debugging)
%
%  .qlook = get_heights structure controlling quick look generation
%   .en = boolean, whether or not to produce a quick look product
%    .out_path = output path of the quick look.  Three forms:
%      1. empty: default path based on the gRadar.out_path, param.records_name,
%         param.season_name, and param.day_seg
%      2. relative path: path based on gRadar.out_path and the contents of
%         .qlook.out_path
%      3. absolute path: uses this specific path for outputs
%   .wf_comb = vector of times of when to combine waveforms
%
% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
% surfTimes = vector of propagation delays to the surface
%
% Author: John Paden
%%
% See also get_heights.m

global g_data;

physical_constants;
surfTimes = [];


  file1 = fullfile(ct_filename_out(param, ...
  param.analysis.out_path,'CSARP_basal_condition'), ...
  sprintf('IceBedCoherenceIndex_%s_%02d_%03d.mat',param.day_seg,param.proc.frm,param.load.blockno));


 if exist(file1,'file')
     success = true;
     return
 end
  

if ~isfield(param.get_heights.qlook,'save_format') || isempty(param.get_heights.qlook.save_format)
  param.get_heights.qlook.save_format = '6';
end

save_format = sprintf('-v%s',param.get_heights.qlook.save_format);

records_fn = ct_filename_support(param,param.records.records_fn,'records');

if ~isfield(param.get_heights,'elev_correction') || isempty(param.get_heights.elev_correction)
  param.get_heights.elev_correction = false;
end

%% Set simple_firdec (boolean, true means decimate in loader for efficiency)
if length(param.get_heights.B_filter) == param.get_heights.decimate_factor ...
    && all(param.get_heights.B_filter == param.get_heights.B_filter(1))
  if ~param.get_heights.elev_correction
    % Most radar headers do not support Elevation1 correction so it must
    % be disabled to allow simple_firdec
    simple_firdec = true;
  elseif any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5'}))
    % FMCW radars have Elevation1 compensation in data loader so they
    % can still have simple_firdec with Elevation1 correction.
    simple_firdec = true;
  else
    simple_firdec = false;
  end
else
  simple_firdec = false;
end

if ~isfield(param.get_heights,'trim_vals') || isempty(param.get_heights.trim_vals)
  param.get_heights.trim_vals = [0 0];
end

if ~isfield(param.get_heights,'coh_noise_method') || isempty(param.get_heights.coh_noise_method)
  param.get_heights.coh_noise_method = 0;
end

if ~isfield(param.get_heights,'coh_noise_arg')
  param.get_heights.coh_noise_arg = [];
end

if ~isfield(param.get_heights,'deconvolution') || isempty(param.get_heights.deconvolution)
  param.get_heights.deconvolution = 0;
end
if ~isfield(param.get_heights,'deconv_enforce_wf_idx') 
  param.get_heights.deconv_enforce_wf_idx = [];
end
if ~isfield(param.get_heights,'deconv_same_twtt_bin') 
  param.get_heights.deconv_same_twtt_bin = [];
end

if ~isfield(param.get_heights,'psd_smooth') || isempty(param.get_heights.psd_smooth)
  param.get_heights.psd_smooth = 0;
end

if ~isfield(param.get_heights,'ft_oversample') || isempty(param.get_heights.ft_oversample)
  param.get_heights.ft_oversample = 1;
end

if ~isfield(param.get_heights,'pulse_rfi') || isempty(param.get_heights.pulse_rfi)
  param.get_heights.pulse_rfi.en = 0;
end

if ~isfield(param.get_heights,'ft_dec') || isempty(param.get_heights.ft_dec)
  param.get_heights.ft_dec = 1;
end

if ~isfield(param.get_heights,'pulse_comp') || isempty(param.get_heights.pulse_comp)
  param.get_heights.pulse_comp = 1;
end

if ~isfield(param.get_heights,'raw_data') || isempty(param.get_heights.raw_data)
  param.get_heights.raw_data = 0;
end

if ~isfield(param.get_heights,'roll_correction') || isempty(param.get_heights.roll_correction)
  param.get_heights.roll_correction = 0;
end

if param.get_heights.roll_correction
  param.get_heights.combine_rx = false;
else
  param.get_heights.combine_rx = true;
end

if ~isfield(param.records,'file_version')
  param.records.file_version = [];
end

if abs(sum(param.get_heights.B_filter)-1) > 1e4*eps
  warning('B_filter weights are not normalized. They must be normalized so normalizing to one now.')
  param.get_heights.B_filter = param.get_heights.B_filter / sum(param.get_heights.B_filter);
end

if ~isfield(param.get_heights,'inc_B_filter') || isempty(param.get_heights.inc_B_filter)
  param.get_heights.inc_B_filter = 1;
end
if abs(sum(param.get_heights.inc_B_filter)-1) > 1e4*eps
  warning('inc_B_filter weights are not normalized. They must be normalized so normalizing to one now.')
  param.get_heights.inc_B_filter = param.get_heights.inc_B_filter / sum(param.get_heights.inc_B_filter);
end

% =====================================================================
% Determine which records to load with load_mcords_data
%
% Load records on either side of the current block, note if at the
% beginning or end of the segment.  Load with minimal presumming.

if simple_firdec
  load_param.load.recs(1) = param.load.recs(1);
  load_param.load.recs(2) = param.load.recs(2);
  records = read_records_aux_files(records_fn,load_param.load.recs);
  old_param_records = records.param_records;
else
  if mod(length(param.get_heights.B_filter)-1,2)
    error('Filter order must be even (e.g. fir1(EVEN_NUMBER,cutoff))');
  end
  filter_order = length(param.get_heights.B_filter) - 1;
  start_buffer = min(filter_order/2,param.load.recs(1)-1);
  load_param.load.recs(1) = param.load.recs(1)-start_buffer;
  load_param.load.recs(2) = param.load.recs(2)+filter_order/2;
  param.load.recs_keep(1) = param.load.recs_keep(1)-start_buffer;
  param.load.recs_keep(2) = param.load.recs_keep(2)+filter_order/2;
  records = read_records_aux_files(records_fn,load_param.load.recs);
  load_param.load.recs(2) = load_param.load.recs(1) + length(records.gps_time) - 1;
  param.load.recs_keep(2) = param.load.recs_keep(1) + length(records.gps_time) - 1;
  old_param_records = records.param_records;
end
old_param_records.gps_source = records.gps_source;

if isfield(param.get_heights,'surface_src') && param.get_heights.surface_src
  %% Get the generic layer data path
  layer_path = fullfile(ct_filename_out(param,'layerData','',0));
  
  %% Load the current frame
  layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm));
  layer = load(layer_fn);
  new_surface_gps_time = layer.GPS_time;
  new_surface = layer.layerData{1}.value{2}.data;
  new_bottom = layer.layerData{2}.value{2}.data;
  %% Get the previous frame if necessary
  if records.gps_time(1) < new_surface_gps_time(1)-1
    layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm-1));
    if exist(layer_fn,'file')
      layer = load(layer_fn);
      new_surface_gps_time = [layer.GPS_time new_surface_gps_time];
      new_surface = [layer.layerData{1}.value{2}.data new_surface];
      new_bottom = [layer.layerData{2}.value{2}.data new_bottom];
    end
  end
  %% Get the next frame if necessary
  if records.gps_time(end) > new_surface_gps_time(end)+1
    layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm+1));
    if exist(layer_fn,'file')
      layer = load(layer_fn);
      new_surface_gps_time = [new_surface_gps_time layer.GPS_time];
      new_surface = [new_surface layer.layerData{1}.value{2}.data];
      new_bottom = [new_bottom layer.layerData{2}.value{2}.data];
    end
  end
  %% Since layer files may have overlapping data, sort it
  [new_surface_gps_time new_surface_idxs] = sort(new_surface_gps_time);
  new_surface = new_surface(new_surface_idxs);
  new_bottom = new_bottom(new_surface_idxs);

  %% Do the interpolation and overwrite the records.surface variable
  new_surface = interp1(new_surface_gps_time,new_surface,records.gps_time,'linear','extrap');
  records.surface = new_surface;
  new_bottom = interp1(new_bottom_gps_time,new_bottom,records.gps_time,'linear','extrap');
  records.bottom = new_bottom;
end

% =====================================================================
% Collect waveform information into one structure
%  (used by load_RADARNAME_data)

if simple_firdec
  param.get_heights.presums = param.get_heights.decimate_factor;
else
  param.get_heights.presums = 1;
end
if strcmpi(param.radar_name,'mcrds')
  [wfs,rec_data_size] = load_mcrds_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.get_heights);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(param.radar_name,{'acords','hfrds','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  [wfs,rec_data_size] = load_mcords_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.get_heights);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(param.radar_name,{'accum'}))
  [wfs,rec_data_size] = load_steppedchirp_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.get_heights);
  load_param.load.rec_data_size = rec_data_size;
elseif strcmpi(param.radar_name,'icards')%add icards radar here--------------QISHI
  param.get_heights.pulse_comp      =      false;
  param.get_heights.ft_dec     =      false;
  param.get_heights.presums    =      1;
  param.get_heights.combine_rx =      false;
  [wfs,rec_data_size] = load_icards_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.get_heights);
  load_param.load.rec_data_size = rec_data_size;
elseif any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5'}))
  [path name] = fileparts(records_fn);
  cdf_fn = fullfile(path, sprintf('%s.nc', name));
  try
    records.settings.nyquist_zone = ncread(cdf_fn,'settings(1).nyquist_zone', ...
      [1 load_param.load.recs(1)],[1 diff(load_param.load.recs([1 end]))+1]);
  end
  try
    records.settings.loopback_mode = ncread(cdf_fn,'settings(1).loopback_mode',[1 load_param.load.recs(1)],[1 1]);
  end
  wfs_idx = find(records.settings.wfs_records <= load_param.load.recs(1),1,'last');
  records.settings.wfs = records.settings.wfs(wfs_idx).wfs;
  wfs = load_fmcw_wfs(records.settings, param, ...
    1:max(old_param_records.records.file.adcs), param.get_heights);
end
% load_wf_cmd=['load_',param.radar_name,'_wfs(records.wfs,param,1:max(old_param_records.file.adcs),param.get_heights)'];
% [wfs,rec_data_size] = eval(load_wf_cmd);
load_param.wfs = wfs;

% =====================================================================
% Collect record file information required for using load_RADAR_NAME_data
%  - Performs mapping between param.rxs and the records file contents
%  - Translates filenames from relative to absolute
%  - Makes filenames have the correct filesep
 
% Create a list of unique adcs required by the imgs list
param.load.adcs = [];
for idx = 1:length(param.load.imgs)
  new_adcs = abs(param.load.imgs{idx}(:,2:2:end)).';
  param.load.adcs = unique(cat(2, new_adcs(:).', param.load.adcs));
end

recs = load_param.load.recs - load_param.load.recs(1) + 1;
if any(strcmpi(param.radar_name,{'hfrds','icards','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  % adc_headers: the actual adc headers that were loaded
  if ~isfield(old_param_records.records.file,'adc_headers') || isempty(old_param_records.records.file.adc_headers)
    old_param_records.records.file.adc_headers = old_param_records.records.file.adcs;
  end
  
  % boards_headers: the boards that the actual adc headers were loaded from
  boards_headers = adc_to_board(param.radar_name,old_param_records.records.file.adcs);
  
  for idx = 1:length(param.load.adcs)
    % adc: the specific ADC we would like to load
    adc = param.load.adcs(idx);
    % adc_idx: the records file index for this adc
    adc_idx = find(old_param_records.records.file.adcs == adc);
    if isempty(adc_idx)
      error('ADC %d not present in records file\n', adc);
    end
    
    % board: the board associated with the ADC we would like to load
    board = adc_to_board(param.radar_name,adc);
    % board_header: the board headers that we will use with this ADC
    board_header = adc_to_board(param.radar_name,old_param_records.records.file.adc_headers(adc_idx));
    % board_idx: the index into the records board list to use
    board_idx = find(board_header == boards_headers);
    
    % Just get the file-information for the records we need
    load_param.load.file_idx{idx} = relative_rec_num_to_file_idx_vector( ...
      load_param.load.recs,records.relative_rec_num{board_idx});
    load_param.load.offset{idx} = records.offset(board_idx,:);
    file_idxs = unique(load_param.load.file_idx{idx});
    
    % Recognize if first record is really from previous file and it is a
    % valid record (i.e. offset does not equal -2^31)
    if sign(load_param.load.offset{idx}(1)) < 0 && load_param.load.offset{idx}(1) ~= -2^31
      file_idxs = [file_idxs(1)-1 file_idxs];
    end
    
    % Just copy the filenames we need
    load_param.load.filenames{idx}(file_idxs) = records.relative_filename{board_idx}(file_idxs);

    % Modify filename according to channel
    for file_idx = 1:length(load_param.load.filenames{idx})
      if ~isequal(old_param_records.records.file.adc_headers,old_param_records.records.file.adcs)
        load_param.load.filenames{idx}{file_idx}(9:10) = sprintf('%02d',board);
      end
    end
    
    filepath = get_segment_file_list(param,adc);
    
    % Convert relative file paths into absolute file paths if required,
    % also corrects filesep (\ and /)
    for file_idx = 1:length(load_param.load.filenames{idx})
      load_param.load.filenames{idx}{file_idx} ...
        = fullfile(filepath,load_param.load.filenames{idx}{file_idx});
    end
  end
  load_param.load.file_version = param.records.file_version;
elseif strcmpi(param.radar_name,'mcrds')
  load_param.load.offset = records.offset;
  load_param.load.file_rec_offset = records.relative_rec_num;
  load_param.load.filenames = records.relative_filename;
  base_dir = ct_filename_data(param,param.vectors.file.base_dir);
  adc_folder_name = param.vectors.file.adc_folder_name;
  load_param.load.filepath = fullfile(base_dir, adc_folder_name);
  load_param.load.wfs = records.settings.wfs;
  load_param.load.wfs_records = records.settings.wfs_records;
elseif any(strcmpi(param.radar_name,{'accum'}))
  load_param.load.offset = records.offset;
  load_param.load.file_rec_offset = records.relative_rec_num;
  load_param.load.filenames = records.relative_filename;
  base_dir = ct_filename_data(param,param.vectors.file.base_dir);
  adc_folder_name = param.vectors.file.adc_folder_name;
  load_param.load.filepath = fullfile(base_dir, adc_folder_name);
  load_param.load.wfs = records.settings.wfs;
  load_param.load.radar_name = param.radar_name;
  load_param.load.season_name = param.season_name;
  load_param.load.day_seg = param.day_seg;
  load_param.load.tmp_path = param.tmp_path;
  load_param.load.file_version = param.records.file_version;
elseif strcmpi(param.radar_name,'acords')
  load_param.load.offset = records.offset;
  load_param.load.file_rec_offset = records.relative_rec_num;
  load_param.load.filenames = records.relative_filename;
  base_dir = ct_filename_data(param,param.vectors.file.base_dir);
  adc_folder_name = param.vectors.file.adc_folder_name;
  load_param.load.filepath = fullfile(base_dir, adc_folder_name);
  load_param.load.wfs = records.settings.wfs;
  load_param.load.wfs_records = records.settings.wfs_records;
elseif any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5'}))
  % Determine which ADC boards are supported and which ones were actually loaded
  if ~isfield(old_param_records.records.file,'adc_headers') || isempty(old_param_records.records.file.adc_headers)
    old_param_records.records.file.adc_headers = old_param_records.records.file.adcs;
  end
  boards = adc_to_board(param.radar_name,old_param_records.records.file.adcs);
  boards_headers = adc_to_board(param.radar_name,old_param_records.records.file.adc_headers);
  
  for idx = 1:length(param.load.adcs)
    adc = param.load.adcs(idx);
    adc_idx = find(old_param_records.records.file.adcs == param.load.adcs(idx));
    if isempty(adc_idx)
      error('ADC %d not present in records file\n', param.load.adcs(idx));
    end
    
    % Just get the file-information for the records we need
    board = adc_to_board(param.radar_name,adc);
    actual_board_idx = find(board == boards);
    board_idx = find(old_param_records.records.file.adc_headers(actual_board_idx) == boards_headers);
    load_param.load.file_idx{idx} = relative_rec_num_to_file_idx_vector( ...
      load_param.load.recs,records.relative_rec_num{board_idx});
    load_param.load.offset{idx} = records.offset(board_idx,:);
    file_idxs = unique(load_param.load.file_idx{idx});
    
    % Recognize if first record is really from previous file and it is a
    % valid record (i.e. offset does not equal -2^31)
    if sign(load_param.load.offset{idx}(1)) < 0 && load_param.load.offset{idx}(1) ~= -2^31
      file_idxs = [file_idxs(1)-1 file_idxs];
    end
    
    % Just copy the filenames we need
    load_param.load.filenames{idx}(file_idxs) = records.relative_filename{board_idx}(file_idxs);
    
    % Modify filename according to channel
    for file_idx = 1:length(load_param.load.filenames{idx})
      if ~isequal(old_param_records.records.file.adc_headers,old_param_records.records.file.adcs)
        load_param.load.filenames{idx}{file_idx}(9:10) = sprintf('%02d',board);
      end
    end
    
    filepath = get_segment_file_list(param,adc);
    
    % Convert relative file paths into absolute file paths if required,
    % also corrects filesep (\ and /)
    for file_idx = 1:length(load_param.load.filenames{idx})
      load_param.load.filenames{idx}{file_idx} ...
        = fullfile(filepath,load_param.load.filenames{idx}{file_idx});
    end
  end
  load_param.load.file_version = param.records.file_version;
else
  error('Radar name %s not supported', param.radar_name);
end

% =====================================================================
% Setup control parameters for load_mcords_data

load_param.load.adcs = param.load.adcs;

load_param.proc.trim_vals           = param.get_heights.trim_vals;
load_param.proc.pulse_comp          = param.get_heights.pulse_comp;
load_param.proc.ft_dec              = param.get_heights.ft_dec;
load_param.proc.ft_wind             = param.get_heights.ft_wind;
load_param.proc.ft_wind_time        = param.get_heights.ft_wind_time;
load_param.proc.presums             = param.get_heights.presums;
load_param.proc.combine_rx          = param.get_heights.combine_rx;
load_param.proc.pulse_rfi           = param.get_heights.pulse_rfi;
load_param.proc.coh_noise_method    = param.get_heights.coh_noise_method;
load_param.proc.coh_noise_arg       = param.get_heights.coh_noise_arg;

load_param.radar = param.radar;
load_param.surface = records.surface;
if strcmpi(param.radar_name,'acords')
  param.gps_source = records.gps_source;
end

% =====================================================================
% Load and process each waveform separately
%
% For each waveform:
% 1. Load receiver data separately (minimal presumming)
% 2. Remove coherent noise (slow-time mean removal)
% 3. Apply roll correction
% 4. Combine receivers
% 5. Apply Elevation1 compensation
% 6. FIR decimate the data
% =====================================================================

for img_idx = 1:length(param.load.imgs)
  % Setup roll correction
  if param.get_heights.roll_correction
    if isempty(param.get_heights.lever_arm_fh)
      error('lever_arm_fh must be defined if roll_correction is enabled');
    end
    
    trajectory_param = struct('gps_source',records.gps_source, ...
      'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
      'tx_weights', [], 'lever_arm_fh', param.get_heights.lever_arm_fh);
    ref = trajectory_with_leverarm(records,trajectory_param);

    lever_arm_fh = param.get_heights.lever_arm_fh;
    % Setup motion compensation (roll removal)
    radar_lever_arm = zeros(3,size(param.load.imgs{img_idx},1));
    for wf_adc_idx = 1:size(param.load.imgs{img_idx},1)
      wf = abs(param.load.imgs{img_idx}(wf_adc_idx,1));
      adc = abs(param.load.imgs{img_idx}(wf_adc_idx,2));
      radar_lever_arm(:,wf_adc_idx) = lever_arm_fh(trajectory_param,wfs(wf).tx_weights,wfs(wf).rx_paths(adc));
    end
  end
  
  % Default values to use
  wf = abs(param.load.imgs{img_idx}(1,1));
  adc = abs(param.load.imgs{img_idx}(1,2));
  lambda_fc = c/wfs(wf).fc;

  %% Compute trajectory using GPS/INS data and the lever arm
  % out_records = motion compensated data
  % records = original record information (not to be used below this section)
  if isempty(param.get_heights.lever_arm_fh)
    out_records = records;
  else
    trajectory_param = struct('gps_source',old_param_records.gps_source, ...
          'season_name',param.season_name,'radar_name',param.radar_name, ...
          'rx_path', wfs(wf).rx_paths(adc), ...
      'tx_weights', wfs(wf).tx_weights, 'lever_arm_fh', param.get_heights.lever_arm_fh);
    for tmp_wf_adc_idx = 2:size(param.load.imgs{1},1)
      tmp_wf = abs(param.load.imgs{1}(tmp_wf_adc_idx,1));
      tmp_adc = abs(param.load.imgs{1}(tmp_wf_adc_idx,2));
      trajectory_param.rx_path(tmp_wf_adc_idx) = wfs(tmp_wf).rx_paths(tmp_adc);
    end
    out_records = trajectory_with_leverarm(records,trajectory_param);
  end
  
  %% Load data into g_data using load_mcords_data
  load_param.load.imgs = param.load.imgs(img_idx);
  % Determine combination times when multiple wf-adc pairs are being loaded
  % to form a single range line
  if size(load_param.load.imgs{1},2) == 2
    load_param.load.wf_adc_comb.en = 0;
  else
    load_param.load.wf_adc_comb.en = 1;
    wf1 = load_param.load.imgs{1}(1,1);
    wf2 = load_param.load.imgs{1}(1,3);
    % t1 = time to switch from wf1 to wf2
    %wf_adc_surface = fir_dec(records.surface, param.get_heights.decimate_factor);
    wf_adc_surface = records.surface;
    t1 = eval(param.get_heights.qlook.wf_adc_comb{1});
    % If the time to switch is longer than wf1 then it gets capped to wf1
    t1(t1 > wfs(wf1).time(end) - wfs(wf1).Tpd) = wfs(wf1).time(end) - wfs(wf1).Tpd;
    load_param.load.wf_adc_comb.Nt = round((wfs(wf2).time(end) - wfs(wf1).time(1)) / wfs(wf1).dt);
    load_param.load.wf_adc_comb.rbins(1,:) = round((t1-wfs(wf1).time(1))/wfs(wf1).dt);
    load_param.load.wf_adc_comb.rbins(2,:) = 1+round((t1-wfs(wf2).time(1))/wfs(wf2).dt);
    % 1:load_param.load.img_comb{1}(1,rline) <- 1:load_param.load.img_comb{1}(1,rline)
    % load_param.load.img_comb{1}(1,rline)+1:end <- load_param.load.img_comb{1}(2,rline):end
  end
  if strcmpi(param.radar_name,'mcords')
    load_mcords_data(load_param);
    g_data = g_data{1};
  elseif any(strcmpi(param.radar_name,{'hfrds','mcords2','mcords3','mcords4','mcords5'}))
    load_mcords2_data(load_param);
    g_data = g_data{1};
  elseif strcmpi(param.radar_name,'mcrds')
    if isfield(records,'adc_phase_corr_deg') && isfield(param.radar,'adc_phase_corr_en') && param.radar.adc_phase_corr_en
      load_param.adc_phase_corr_deg = records.adc_phase_corr_deg;
    else
      load_param.adc_phase_corr_deg = zeros(length(load_param.surface),max(records.param_records.records.file.adcs));
    end
    load_mcrds_data(load_param);
    g_data = g_data{1};
  elseif strcmpi(param.radar_name,'acords')
    load_param.load.file_version = param.records.file_version;
    load_acords_data(load_param);
    g_data = g_data{1};
  elseif strcmpi(param.radar_name,'icards')%add icards----qishi
    load_icards_data(load_param,param);
    g_data = g_data{1};
  elseif strcmpi(param.radar_name,'accum')
    %load_param.proc.elev_correction = param.get_heights.elev_correction;
    load_param.proc.elev_correction = 0;
    load_param.radar_name = param.radar_name;
    load_param.season_name = param.season_name;
    load_param.day_seg = param.day_seg;
    load_param.out_path = param.out_path;
    [g_data,wfs.time] = load_accum_data(load_param,out_records);
    if param.get_heights.ft_oversample ~= 1
      dt = (wfs.time(2)-wfs.time(1))/param.get_heights.ft_oversample;
      Nt = size(g_data,1)*param.get_heights.ft_oversample;
      g_data = interpft(g_data,Nt);
      wfs.time = wfs.time(1) + dt*(0:Nt-1).';
    end
  elseif strcmpi(param.radar_name,'accum2')
    load_accum2_data(load_param);
    g_data = g_data{1};
  elseif any(strcmpi(param.radar_name,{'kuband','snow','kuband2','snow2','kuband3','snow3','kaband3','snow5'}))
    load_param.proc.elev_correction = param.get_heights.elev_correction;
    load_param.proc.deconvolution = param.get_heights.deconvolution;
    load_param.proc.deconv_enforce_wf_idx = param.get_heights.deconv_enforce_wf_idx;
    load_param.proc.deconv_same_twtt_bin = param.get_heights.deconv_same_twtt_bin;
    load_param.proc.psd_smooth = param.get_heights.psd_smooth;
    load_param.radar_name = param.radar_name;
    load_param.season_name = param.season_name;
    load_param.support_path = param.support_path;
    load_param.day_seg = param.day_seg;
    load_param.load.tmp_path = param.tmp_path;
    load_param.out_path = param.out_path;
    [img_time,img_valid_rng,img_deconv_filter_idx,img_Mt] = load_fmcw_data(load_param,out_records);
    g_data = g_data{1};
    valid_rng = img_valid_rng{1};
    deconv_filter_idx = img_deconv_filter_idx{1};
    % Get heights only loads one image at a time, so img_time{1}
    for wf = 1:length(wfs)
      wfs(wf).time = img_time{1};
    end
    if param.get_heights.ft_oversample ~= 1
      dt = (wfs(wf).time(2)-wfs(wf).time(1))/param.get_heights.ft_oversample;
      Nt = size(g_data,1)*param.get_heights.ft_oversample;
      g_data = ifft(fft(g_data),Nt);
      wfs(wf).time = wfs(wf).time(1) + dt*(0:Nt-1).';
      valid_rng = valid_rng*param.get_heights.ft_oversample;
      for rline=1:size(g_data,2)
        g_data([1:valid_rng(1,rline)-1,valid_rng(2,rline)+1:end],rline) = 0;
      end
    end
  end
  if ~exist('deconv_filter_idx','var')
    deconv_filter_idx = NaN(1,size(g_data,2));
  end
  
  if 0
    % Deconvolution Test Code
    keyboard
    dd = load('/mnt/products/ref_Tpd_3us_adc_1.mat');
    %                 ee = fft(dd.ref(1:20:end),length(g_data(:,1,1)));
    dd.ref = ifft(dd.fref);
    ee = interpft(fft(dd.ref(1:20:end)),length(g_data(:,80,1)));
    plot(lp(ee))
    ff = g_data(:,80,1);
    
    figure(1); clf;
    plot(lp(ff));
    hold on
    plot(lp(ifft(fft(ff) .* ee)),'r')
    hold off;
    
    ff = g_data(:,:,1);
    
    figure(1); clf;
    imagesc(lp(ff));
    figure(2); clf;
    imagesc(lp(ifft(fft(ff) .* repmat(ee,[1 size(ff,2)]))))
  end
  
  %% Remove overlap data
  recs_keep = 1+param.load.recs_keep(1)-load_param.load.recs(1) ...
    : length(out_records.lat)+param.load.recs_keep(end)-load_param.load.recs(end);
  out_records.lat = out_records.lat(recs_keep);
  out_records.lon = out_records.lon(recs_keep);
  out_records.elev = out_records.elev(recs_keep);
  out_records.roll = out_records.roll(recs_keep);
  out_records.pitch = out_records.pitch(recs_keep);
  out_records.heading = out_records.heading(recs_keep);
  out_records.gps_time = out_records.gps_time(recs_keep);
  out_records.surface = out_records.surface(recs_keep);
  if simple_firdec
    recs_keep = 1+(param.load.recs_keep(1)-load_param.load.recs(1))/param.get_heights.decimate_factor ...
      : size(g_data,2)+(param.load.recs_keep(end)-load_param.load.recs(end))/param.get_heights.decimate_factor;
  end
  g_data = g_data(:,recs_keep,:);
  deconv_filter_idx = deconv_filter_idx(:,recs_keep);
  
  %% Remove coherent noise
  if param.get_heights.coh_noise_method && ~any(strcmpi(param.radar_name,{'kuband','snow','kuband2','snow2','kuband3','kaband3','snow3','snow5'}))
    
    if param.get_heights.coh_noise_method == 3 && isempty(param.get_heights.coh_noise_arg)
      param.get_heights.coh_noise_arg = 255;
    end
    
    % Remove only the DC Doppler component
    for wf_adc_idx = 1:size(g_data,3)
      if param.get_heights.coh_noise_method == 1
        g_data(:,:,wf_adc_idx) = bsxfun(@minus, g_data(:,:,wf_adc_idx), ...
          mean(g_data(:,:,wf_adc_idx),2));
      elseif param.get_heights.coh_noise_method == 3
        g_data(:,:,wf_adc_idx) = bsxfun(@minus, g_data(:,:,wf_adc_idx), ...
          fir_dec(g_data(:,:,wf_adc_idx),hanning(param.get_heights.coh_noise_arg).'/(param.get_heights.coh_noise_arg/2+0.5),1));
      else
        error('param.get_heights.coh_noise_method %d not supported.',param.get_heights.coh_noise_method);
      end
    end
    
  end

  %% Roll compensation
  if param.get_heights.roll_correction
    % Apply roll-only motion compensation
    for wf_adc_idx = 1:size(g_data,3)
      wf = abs(param.load.imgs{img_idx}(wf_adc_idx,1));
      adc = abs(param.load.imgs{img_idx}(wf_adc_idx,2));
      rx = wfs(wf).rx_paths(adc);
      for rline = 1:size(g_data,2)
        drange = radar_lever_arm(2,wf_adc_idx) * -tan(out_records.roll(rline));
        dphase = drange * 2 * 2*pi/lambda_fc;
        g_data(:,rline,wf_adc_idx) = g_data(:,rline,wf_adc_idx) * exp(1j*dphase);
      end
    end
    g_data = mean(g_data,3);
  end
  
  %% Elevation1 compensation
  if param.get_heights.elev_correction && ~simple_firdec
    % Remove Elevation1 variations (just a phase shift, not a time shift)
    %  - With simple_firdec there is no point in Elevation1 compensation
    %    because the along-track averaging has already been done
    drange = out_records.elev-mean(out_records.elev);
    dphase = drange * 2 * 2*pi/lambda_fc;
    for rline = 1:size(g_data,2)
      g_data(:,rline) = g_data(:,rline) .* exp(1j*dphase(rline));
    end
  end

  %% FIR Decimate
  if simple_firdec
%     if img_idx == 1
      out_records.gps_time = fir_dec(out_records.gps_time, param.get_heights.decimate_factor);
      out_records.lat = fir_dec(out_records.lat, param.get_heights.decimate_factor);
      out_records.lon = fir_dec(out_records.lon, param.get_heights.decimate_factor);
      out_records.elev = fir_dec(out_records.elev, param.get_heights.decimate_factor);
      out_records.roll = fir_dec(out_records.roll, param.get_heights.decimate_factor);
      out_records.pitch = fir_dec(out_records.pitch, param.get_heights.decimate_factor);
      out_records.heading = fir_dec(out_records.heading, param.get_heights.decimate_factor);
%     end
    
  else
    % Check for edge conditions that caused not enough data to be loaded
    % in the case where a segment starts and stops.
    % If not enough data was loaded, modify the filter coefficient
    % normalization so that it handles this
    Nidxs = floor((param.load.recs(2)-param.load.recs(1)+1) / param.get_heights.decimate_factor);
    rline0 = 1 + start_buffer;
    g_data = fir_dec(g_data, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);

%     if img_idx == 1
      out_records.gps_time = fir_dec(out_records.gps_time, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.lat = fir_dec(out_records.lat, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.lon = fir_dec(out_records.lon, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.elev = fir_dec(out_records.elev, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.roll = fir_dec(out_records.roll, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.pitch = fir_dec(out_records.pitch, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
      out_records.heading = fir_dec(out_records.heading, param.get_heights.B_filter, ...
        param.get_heights.decimate_factor, rline0, Nidxs);
%     end
  end
  
  if 0
    % Enable this if-statement only for debugging
    figure(img_idx); clf;
    imagesc([1 size(g_data,2)],wfs(img_idx).time, ...
      lp(g_data));
    keyboard
  end
  
  %% Apply incoherent averaging with decimation
  if size(param.get_heights.inc_B_filter,2) >= 1
    data_incoh = [];
    for adc_idx = 1:size(g_data,3)
      data_incoh(:,:,adc_idx) = fir_dec(fir_dec(abs(g_data(:,:,adc_idx)).^2,param.get_heights.inc_B_filter,1), param.get_heights.inc_ave);
    end
  end
  
  if 0
    % Undo Elevation1 phase correction
    for rline = 1:size(g_data,2)
      drange = out_records.elev-mean(out_records.elev);
      dphase = drange * 2 * 2*pi/lambda_fc;
      g_data(:,rline) = g_data(:,rline) .* exp(-1j*dphase(rline));
    end
    
    % Leading edge detector
    data_lp = lp(data_incoh);
    [maxVal maxIdx] = max(data_lp);
    sizeBins = zeros(size(maxVal));
    for rline = 1:size(data_lp,2)
      surfBins(rline) = find(data_lp(:,rline) > maxVal(rline)-10,1);
      [tmp maxOffset] = max(data_lp(surfBins(rline)+(0:4),rline));
      surfBins2(rline) = surfBins(rline) + maxOffset-1;
    end
    surfBins3 = medfilt1(surfBins2,5);
    surfBins4 = reshape(repmat(surfBins3,[param.get_heights.inc_ave 1]),[1 size(g_data,2)]);
    
    % Phase and Elevation1 comparison (used for finding UTC time offset
    % errors)
    physical_constants;
    
    figure(1); clf;
    plot(out_records.elev);
    hold on;
    plot(-unwrap(angle(g_data(179,:))) * lambda_fc/(2*pi)/2 + 2463-0.1,'r');
    hold off;
    keyboard
  end
  
  % Save quick look results
  
  if param.get_heights.inc_ave >= 1
    Data = data_incoh;
  else
    Data = g_data;
  end
  
  Time = wfs(wf).time;
  
  GPS_time = fir_dec(out_records.gps_time,param.get_heights.inc_ave);
  Latitude = fir_dec(out_records.lat,param.get_heights.inc_ave);
  Longitude = fir_dec(out_records.lon,param.get_heights.inc_ave);
  Elevation = fir_dec(out_records.elev,param.get_heights.inc_ave);
  Roll = fir_dec(out_records.roll,param.get_heights.inc_ave);
  Pitch = fir_dec(out_records.pitch,param.get_heights.inc_ave);
  Heading = fir_dec(out_records.heading,param.get_heights.inc_ave);
  deconv_filter_idx = fir_dec(deconv_filter_idx,param.get_heights.inc_ave);
  
  if ~isnan(out_records.gps_time(1))
    fn = fullfile(ct_filename_out(param, ...
      param.get_heights.qlook.out_path, 'CSARP_qlook'), ...
      sprintf('ql_data_%03d_01_01',param.proc.frm), sprintf('%s_img_%02d.mat', ...
      datestr(epoch_to_datenum(out_records.gps_time(1)), 'yyyymmdd_HHMMSS'), ...
      img_idx));
  else % added for loopback test data in lab with gps data
    fn = fullfile(ct_filename_out(param, ...
      param.get_heights.qlook.out_path, 'CSARP_qlook'), ...
      sprintf('ql_data_%03d_01_01',param.proc.frm), sprintf('%s_img_%02d.mat', ...
      datestr(now, 'yyyymmdd_HHMMSS'), ...
      img_idx));
  end
  [path name ext] = fileparts(fn);
  if ~exist(path,'dir')
    mkdir(path);
  end
  param_records = old_param_records;
  param_get_heights = param;
  custom = [];
  if any(~isnan(deconv_filter_idx))
    custom.deconv_filter_idx = deconv_filter_idx;
  end
  clear deconv_filter_idx;
  
  
  
 %Trim time<0
      first_idx = find(Time <= 0,1,'last');
      if ~isempty(first_idx)
        Time = Time(first_idx:end);
        Data = Data(first_idx:end,:,:);
      end
  
      
      %%  
      
layer = load(fullfile(ct_filename_out(param,'CSARP_post/layerData',''), ...
      sprintf('Data_%s_%03d.mat',param.day_seg,param.proc.frm)));
distance=geodetic_to_along_track( Latitude, Longitude, Elevation); %Horizontal Distance
sf = interp1(layer.GPS_time,layer.layerData{1}.value{2}.data, GPS_time,'linear','extrap'); %Two way prop time to surface
bttime = interp1(layer.GPS_time, layer.layerData{2}.value{2}.data, GPS_time,'linear','extrap');  %Two way prop. time to Bottom
bttime(bttime> Time(end) | bttime< Time(1)) = NaN; % remove pick errors exist
 

power=lp(Data);
 bttime_n=NaN*ones(1,size(Data,2));  
for rline=1:size(Data,2)
   if isnan(bttime(rline))
       continue;
   end

    bt_idx_m = find(Time>bttime(rline),1,'first'); 
    if bt_idx_m-15>0 && bt_idx_m+15<size(Data,1)
        [~,bt_idx]=max(power(bt_idx_m-15:bt_idx_m+15,rline));
        dt=Time(2)-Time(1);
        bttime_n(rline)=(bttime(rline)-15*dt)+(bt_idx-1)*dt;
    if bttime_n(rline)<sf(rline)
        bttime_n(rline)=sf(rline);
    end
   end
end

 debug_flag=0;
  if debug_flag==1
      figure(1); imagesc(distance, Time,lp( Data));
      hold on; plot(distance,sf,'--'); plot(distance,bttime_n,'--');
  end;
 %% 
p=4.99;  %Pre windowed radar pulse halfwidth in air    
slopeval=zeros(1,size( Data,1));

bttime_temp=bttime;
sf_temp=sf;
for rline=1:size(bttime,2)
if isnan(bttime_temp(rline))
    sf_temp(rline)=NaN;
    
    sf(rline)=NaN;
end
if isnan(sf_temp(rline))
   bttime_temp(rline)=NaN; 
   bttime_n(rline)=NaN;
end
if isinf(sf_temp(rline))
    sf_temp(rline)=NaN;
    sf(rline)=NaN;
end
if isinf(bttime_temp(rline))
    bttime_temp(rline)=NaN;
    bttime_n(rline)=NaN;
end 
end
sf(sf==Inf)=NaN;
MeanDepth_og=(nanmean(bttime_temp)-nanmean(sf_temp))*c/(2*sqrt(er_ice));  %Original Mean Depth
MeanDepth=(nanmean(bttime_n)-nanmean(sf))*c/(2*sqrt(er_ice));    %Mean Depth

if (MeanDepth_og~=0 && (~isnan(MeanDepth_og))) && MeanDepth_og>150
    Nx_int_dist=sqrt(MeanDepth*p/sqrt(er_ice));                %Coherent Integration Distance
    Nx_int = floor(2*Nx_int_dist/(distance(10)-distance(9)));
    if Nx_int==0
        continue;
    end
    %% Correction for Elevation1
    sf_new = zeros(size(sf));
    bt_new = zeros(size(bttime_n));
    Elevation1_new = zeros(size( Elevation));
    sf_elev_new = zeros(size( Elevation));
    bt_elev_new = zeros(size( Elevation));
    
    % 1)Remove data before zero time
    negative_bins =  Time < 0;
    Time_new =  Time(~negative_bins);
    Data_new =  Data(~negative_bins,:);
    
    % 2)Create Elevation1 axis to interpolate to
    [max_elev,max_elev_idx] = max( Elevation);
    min_elev = min( Elevation - sf*c/2 - ( Time_new(end)-sf)*c/2/sqrt(er_ice));
    dt =  Time(2)- Time(1);
    dr = dt * c/2 / sqrt(er_ice);
    dt_air = dr/(c/2);
    elev_axis = max_elev:-dr:min_elev;
    new_time = zeros(length(elev_axis),length( Elevation));
    
    % 3)Zero pad data to create space for interpolated data
    zero_pad_len = length(elev_axis) - length( Time_new);
    Data_new = cat(1, Data_new,zeros(zero_pad_len,size( Data_new,2)));
    
    % 4)Determine the corrections to apply to Elevation1 and layers
    dRange = max_elev -  Elevation;
    dBins = round(dRange / (c/2) / dt);
    dtime = dRange/(c/2);
    
    for rline = 1:size( Data_new,2)
        % Determine Elevation1 bins before surface
        sf_elev =  Elevation(rline) - sf(rline) * c/2;
        time0 = -(max_elev -  Elevation(rline))/(c/2);
        last_air_idx = find(elev_axis > sf_elev,1,'last');
        new_time_tmp = (time0 + dt_air*(0:last_air_idx-1)).';
        if last_air_idx < length(elev_axis)
            % Determine Elevation1 bins after surface
            dt_ice = dr/(c/2/sqrt(er_ice));
            first_ice_idx = last_air_idx + 1;
            time0 = sf(rline) + (sf_elev - elev_axis(first_ice_idx))/(c/2/sqrt(er_ice));
            new_time(:,rline) = cat(1,new_time_tmp, (time0 + dt_ice*(0:length(elev_axis)-length(new_time_tmp)-1)).');
        end
        Data_new(:,rline) = interp1( Time_new,  Data_new(1:length( Time_new),rline), new_time(:,rline), 'linear',0);
        Elevation1_new(rline) =  Elevation(rline) + dRange(rline);
        sf_new(rline) = sf(rline) + dtime(rline);
        bt_new(rline) = bttime_n(rline) + dtime(rline);
        sf_elev_new(rline) =  Elevation1_new(rline) - sf_new(rline)*c/2;
        bt_elev_new(rline) =  sf_elev_new(rline) - (bt_new(rline)-sf_new(rline))*c/2/sqrt(er_ice);
    end
    
    if debug_flag==1
        fh = figure(2);
        figure(fh);imagesc([],elev_axis,10*log10(abs( Data_new).^2));title('Elevation1 Correction Data')
        ax = gca;
        ax.YDir = 'normal';
        hold on;plot( sf_elev_new,'--');plot( bt_elev_new,'--');
        bt_slope = diff( bt_elev_new)./diff(distance);
    end
    
    %% Ice Slope Correction
    
    Nx0 = size( Data_new,2);
    Nx = floor(Nx0/Nx_int);
    Nx_mod = mod(Nx0,Nx_int);
    if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
    end
    
    for rline = 1:Nx
        
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        if isnan(rline)
            continue;
        end
        if isinf(mean( bt_elev_new(idx1:idx2))) || isnan(mean( bt_elev_new(idx1:idx2))) % If No Ice Bottom, skip
            continue;
        end
        
        p = polyfit(distance(idx1:idx2), bt_elev_new(idx1:idx2),1);    % Polygonal fitting
        slopeval(idx1:idx2) = polyval(p,distance(idx1:idx2));           % Polygonal values after fitting
        base = distance(idx2)-distance(idx1);
        perpendicular = slopeval(idx2)-slopeval(idx1);
        hypotenuse = sqrt(base^2+perpendicular^2);
        angle1(rline)=asin(perpendicular/hypotenuse)*180/pi;
        slopeerror = slopeval(idx1:idx2)-slopeval(idx1);         % Error betn original and fitting line
        dtime = 2*slopeerror/c/sqrt(er_ice);
        
        for idx = idx1:idx2
            Data_new(:,idx) = interp1(elev_axis,  Data_new(:,idx), elev_axis + slopeerror(idx-idx1 +1), 'linear',0);
            Elevation1_new(idx) =  Elevation1_new(idx) - slopeerror(idx-idx1 +1);
            sf_elev_new(idx) =  sf_elev_new(idx) - slopeerror(idx-idx1 +1);
            bt_elev_new(idx) =  bt_elev_new(idx) - slopeerror(idx-idx1 +1);
            sf_new(idx) = sf_new(idx) + dtime(idx-idx1 +1);
            bt_new(idx) = bt_new(idx) + dtime(idx-idx1 +1);
        end
        
    end
    
    if debug_flag==1
        fh = figure(3);
        figure(fh);imagesc([],elev_axis,10*log10(abs( Data_new).^2));title('Slope Correction Data')
        ax = gca;
        ax.YDir = 'normal';
        hold on;plot( sf_elev_new,'--');plot( bt_elev_new,'--');
    end
    
    %%
    % Truncate data around ice bottom within bt.range_bins
    bt_range_bins =param.analysis.bt.range_bins;
    bt.peakval = NaN*ones(1,size( Data_new,2));
    bt.peakidx = NaN*ones(1,size( Data_new,2));
    bt.waveform = NaN*ones(length(bt_range_bins),size( Data_new,2));
    bt.inc_wf_ave=zeros(length(bt_range_bins),Nx);
    for rline = 1:size( Data_new,2)
        if ~isnan( bt_elev_new(rline)) & ~isinf( bt_elev_new(rline))
            % bt.idx(rline) = find(elev_axis<= bt_elev_new(rline),1,'first');
            bt.peakidx(rline) = round(interp1(elev_axis,[1:length(elev_axis)], bt_elev_new(rline)));
            bt.peakval(rline) =  Data_new(bt.peakidx(rline),rline);
            first_idx = bt.peakidx(rline) + bt_range_bins(1);
            last_idx = bt.peakidx(rline) + bt_range_bins(end);
            if first_idx < 1 | last_idx>size( Data_new,1)
                bt.peakidx(rline) = NaN;
                bt.peakval(rline) = NaN;
                continue
            end
            lower=bt.peakidx(rline)+bt_range_bins(1);
            upper=bt.peakidx(rline)+bt_range_bins(end);
            if upper >size( Data_new,1)
                upper=size( Data_new,1);
            end
            bt.waveform(1:51+upper-bt.peakidx(rline),rline) =  Data_new(lower:upper, rline);
        else
            continue
        end
    end
    
    if debug_flag==1
        figure(4);imagesc(lp( Data_new));
        hold on;plot(bt.peakidx,'--');
        figure(5);plot(lp(bt.waveform));
        figure(6);imagesc(lp(bt.waveform));
    end
    
    
    %% Coherene Index Calculation with Slope error Correction
    
    bt.coherenceindex.Latitude_mean=NaN*ones(1,Nx);
    bt.coherenceindex.Longitude_mean=NaN*ones(1,Nx);
    bt.coherenceindex.GPS_time_ave=NaN*ones(1,Nx);
    square_int = zeros(size(  Data_new,1),Nx); % incoherent integration, take square of abs first, then sum
    square_int_dB=zeros(size(  Data_new,1),Nx);
    int_square = zeros(size( Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
    coh_index_slopecorr = NaN*ones(1,Nx);
    Abruptiveindex=NaN*ones(1,Nx);
    risingedge=NaN*ones(1,Nx);
    fallingedge=NaN*ones(1,Nx);
    bt.avepeakindex=NaN*ones(1,Nx);
    Padj=NaN*ones(1,Nx);
    for rline =1:Nx
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        
        
        square_int(:,rline) = mean(abs( Data_new(:,idx1:idx2)).^2,2);
        int_square(:,rline) = abs(mean( Data_new(:,idx1:idx2),2)).^2;
        square_int_dB(:,rline) = 10*log10(square_int(:,rline));
        
        meanbt=nanmean( bt_elev_new(idx1:idx2));
        if meanbt==0 | isnan(meanbt)
            continue;                %skip if no ice bottom
        end
        
        meansf=nanmean(sf_elev_new(idx1:idx2));   %If bottm=surface then skip
        if meanbt==meansf
            continue;
        end
         depth=meansf-meanbt;
         if depth<150
             continue;
         end
        %bt_idx_m = find(elev_axis<=meanbt,1,'first');
        bt_idx_m = round(interp1(elev_axis,[1:length(elev_axis)],meanbt));
        
        if isempty(bt_idx_m)
            continue;      %skip if no ice bottom
        end
        
        b1=bt_idx_m+param.analysis.bt.peak_bins(1);               % Peak search index peak bins
        b2=bt_idx_m+param.analysis.bt.peak_bins(2);
        if b2>size( Data_new,1)
            b2=size( Data_new,1);
        end
        
        [bt_val,bt_idx] = max(square_int(b1:b2,rline));   %Peak Index and Value
        bt_idx = bt_idx_m+param.analysis.bt.peak_bins(1)+bt_idx-1;
        bt_pwr = 10*log10(bt_val);
        noise_bin1=bt_idx+270;
        noise_bin2=bt_idx+270+param.analysis.bt.noise_bins;
        if noise_bin1>size( Data_new,1) || noise_bin2>size( Data_new,1)
            noise_bin1=size( Data_new,1)-param.analysis.bt.noise_bins;
            noise_bin2=size( Data_new,1);
        end
        
        noise = mean(square_int_dB(noise_bin1:noise_bin2,rline));
        if bt_pwr-noise<3
            
            continue;
        end
        %noise = mean(square_int_dB(noise_bin1:noise_bin2,rline));
        risingedge(rline) = bt_idx-1;
        SNR = bt_pwr-noise;
        
        % find risingedge(rline): risingedge and fallingedge(rline):fallingedge for integration in range
        while square_int_dB(risingedge(rline),rline)-noise > 0.05*SNR  && risingedge(rline)>2
            risingedge(rline) = risingedge(rline) - 1;
        end
        if bt_idx-risingedge(rline)>param.analysis.bt.max_index_range_rising
            risingedge(rline)=bt_idx-param.analysis.bt.max_index_range_rising;  %Max 50 bins away for risingedge
        end
        
        if bt_idx-risingedge(rline)<param.analysis.bt.guard_bins
            risingedge(rline)=bt_idx-param.analysis.bt.guard_bins;  %Guard band of 3
        end
        
        [r,tmp_idx]=(min(square_int_dB(risingedge(rline):bt_idx-param.analysis.bt.guard_bins,rline)));
        risingedge=risingedge+tmp_idx-1;
        
        fallingedge(rline) = bt_idx+1;
        if fallingedge(rline)>size( Data_new,1)
            fallingedge(rline)=size( Data_new,1);
        end
        while square_int_dB(fallingedge(rline),rline)-noise > 0.05*SNR && fallingedge(rline)<size( Data,1)
            fallingedge(rline) = fallingedge(rline) + 1;
        end
        if fallingedge(rline)-bt_idx>param.analysis.bt.max_index_range_falling
            fallingedge(rline)=bt_idx+param.analysis.bt.max_index_range_falling;  %Max 50 bins away for falling edge
        end
        if fallingedge(rline)-bt_idx<param.analysis.bt.guard_bins
            fallingedge(rline)=bt_idx+param.analysis.bt.guard_bins;
        end
        if fallingedge(rline)>size( Data_new,1)
            fallingedge(rline)=size( Data_new,1);
        end
        coh_index_slopecorr(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline)); %Coherence Index
        
        Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));  %Power summed from risingedge(rline) to fallingedge(rline)
        Abruptiveindex(rline)=bt_val/Imeanx;    %Abruptive Index
        
        bt.avepeakindex(rline)=bt_idx;           %peakindex
       
        depth=meansf-meanbt;
        B=2.3*3000/(depth+2000);
        Padj(rline)=(lp(Imeanx*depth.^2)+B*depth/100);
        
        
        bt.coherenceindex.Latitude_mean(rline)=mean( Latitude(idx1:idx2));  %Mean Latitude
        bt.coherenceindex.Longitude_mean(rline)=mean( Longitude(idx1:idx2));
        bt.coherenceindex.GPS_time_ave(rline)=mean( GPS_time(idx1:idx2));
        
        
        if bt_idx-50<=0 | bt_idx+100>size(square_int_dB,1)
            continue;
        end
        bt.inc_wf_ave(:,rline) = square_int_dB(bt_idx+bt_range_bins);  %Average Waveform
    end
   
    if debug_flag==1
        hold on; figure(8); hold on; plot(coh_index_slopecorr,'g','Displayname','SlopeCorrected'); title('Coherence Index')
        legend('show');
        figure(7);plot(coh_index_slopecorr),title('Coherence Index after Slope Correction')
    end
    
    %%
    
    
    bt.coherenceindex.GPS_time= GPS_time;
    bt.coherenceindex.Latitude= Latitude;
    bt.coherenceindex.Longitude= Longitude;
    bt.coherenceindex.Elevation_new= Elevation1_new;
    bt.coherenceindex.Surface_new= sf_elev_new;
    bt.coherenceindex.SurfaceTime_new=sf_new;
    bt.coherenceindex.Bottom_new= bt_elev_new;
    bt.coherenceindex.BottomTime_new=bt_new;
    bt.coherenceindex.Time_new= Time_new';
    bt.coherenceindex.frm_id={};
    bt.coherenceindex.value=coh_index_slopecorr;
    bt.coherenceindex.Abruptiveindex=Abruptiveindex;
    bt.coherenceindex.AdjustedIntensity=Padj;
end  
if MeanDepth_og==0 || isnan(MeanDepth_og) || MeanDepth_og<150
      bt.coherenceindex.GPS_time= [];
    bt.coherenceindex.Latitude= [];
    bt.coherenceindex.Longitude= [];
    bt.coherenceindex.Elevation_new= [];
    bt.coherenceindex.Surface_new= [];
    bt.coherenceindex.SurfaceTime_new=[];
    bt.coherenceindex.Bottom_new= [];
    bt.coherenceindex.BottomTime_new=[];
    bt.coherenceindex.Time_new= [];
    bt.coherenceindex.frm_id={};
    bt.coherenceindex.value=[];
    bt.coherenceindex.Abruptiveindex=[];
    bt.coherenceindex.AdjustedIntensity=[];
   bt.coherenceindex.Latitude_mean=[];
    bt.coherenceindex.Longitude_mean=[];
   bt.coherenceindex.GPS_mean=[];
   bt.peakval = [];
   bt.peakidx = [];
   bt.waveform = [];
   bt.inc_wf_ave= [];
   bt.avepeakindex= [];
end
    %%
    out_fn = fullfile(ct_filename_out(param, ...
        param.analysis.out_path,'CSARP_basal_condition'), ...
        sprintf('IceBedCoherenceIndex_%s_%03d_%02d.mat',param.day_seg,param.proc.frm,param.load.blockno));
    
    [out_fn_dir] = fileparts(out_fn);
    if ~exist(out_fn_dir,'dir')
        mkdir(out_fn_dir);
    end
    
    param_analysis = param;
    fprintf('  Saving outputs %s\n', out_fn);
    save(save_format, out_fn, 'bt','param_analysis');

 
      
      
  
end

success = true;

return;
