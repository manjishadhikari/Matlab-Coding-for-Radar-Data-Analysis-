function IceBedCoherenceIndex_tracker(params,param_override)

% IceBedCoherenceIndex_tracker(param,param_override)
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
% See run_IceBedCoherenceIndex_tracker.m for how to run this function directly.
%
% Authors: Jilu Li
%
% See also: master.m, IceBedCoherenceIndex_tracker_task.m

% =====================================================================
% General Setup
% =====================================================================


if ~isstruct(params)
  % Functional form
  params();
end

for param_idx = 1:length(params)
  param = params(param_idx);
  param = merge_structs(param, param_override);
  
  if isfield(param.cmd,'generic') && param.cmd.generic
    if param.analysis.IceBedCoherenceIndex.en
      dbstack_info = dbstack;
      fprintf('=====================================================================\n');
      fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
      fprintf('=====================================================================\n');
      
      if ~isfield(param.analysis,'out_path')
        param.analysis.out_path = '';
      end
      
      if ~isfield(param.analysis,'in_path')
        param.analysis.in_path = '';
      end
      
      frm_fn = param.analysis.in_path;
      
      % =====================================================================
      % Setup processing
      % =====================================================================
      
      % Get WGS84 ellipsoid parameters
      physical_constants;
      
      % Load frames file
      
      load(ct_filename_support(param,param.records.frames_fn,'frames'));
%       layer = load(ct_filename_out(param,param.records.frames_fn,'CSARP_layerData'));
      
      if isempty(param.cmd.frms)
        warning('All frames are processed, setting param.cmd.frms to do all frames.');
        param.cmd.frms = (1:length(frames.frame_idxs)); % All frames
      end
      
      global g_data;
      g_data = [];
      
  
      
      analysis_out_path = ct_filename_out(param, param.analysis.out_path, 'CSARP_basal_condition');
      
      if ~isfield(param.sched,'rerun_only') || isempty(param.sched.rerun_only)
        param.sched.rerun_only = false;
      end
      
      seg_file = fullfile(ct_filename_out(param,param.analysis.out_path, 'CSARP_basal_condition'), ...
            sprintf('IceBedCoherenceIndex_%s.mat',param.day_seg));
        seg_file= regexprep(seg_file,['/',param.day_seg],'');
        existdir=0;
        if exist(seg_file,'file')
            existdir = 1;
        end
        
        if existdir==0
            
            
            
            
            
            % Cleanup files in folders
            if ~param.sched.rerun_only
                if exist(analysis_out_path,'dir')
                    warning('Removing IceBedCoherenceIndex*.mat files in path %s', analysis_out_path);
                    delete(sprintf('%s/IceBedCoherenceIndex*.mat',analysis_out_path));
                end
            end
            
            % =====================================================================
            % Setup static inputs for IceBedCoherenceIndex_tracker_task
            % =====================================================================
            task_param = param;
            
            % =====================================================================
            % Setup the scheduler
            % =====================================================================
            
            if strcmpi(param.sched.type,'custom_torque')
                global ctrl; % Make this global for convenience in debugging
                ctrl = torque_new_batch(param);
                fprintf('Torque batch: %s\n', ctrl.batch_dir);
                torque_compile('IceBedCoherenceIndex_tracker_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
                
            elseif ~strcmpi(param.sched.type,'no scheduler')
                fd = [get_filenames(param.path,'','','.m',struct('recursive',1)); ...
                    get_filenames(param.path,'','','.mexa64',struct('recursive',1))];
                fd_override = [get_filenames(param.path_override,'','','.m',struct('recursive',1)); ...
                    get_filenames(param.path_override,'','','.mexa64',struct('recursive',1))];
                
                fd = merge_filelists(fd, fd_override);
                
                % Remove SVN files from path
                non_svn_path = zeros(size(fd));
                for fn_idx = 1:length(fd)
                    if isempty(strfind(fd{fn_idx},'.svn'))
                        non_svn_path(fn_idx) = 1;
                    end
                end
                fd = fd(non_svn_path);
                
                % Initialize submission ctrl structure
                global ctrl; % Make this global for convenience in debugging
                ctrl = [];
                ctrl.cmd = 'init';
                ctrl.sched = param.sched;
                ctrl.fd = fd;
                ctrl = create_task(ctrl);
                
                % Prepare submission ctrl structure for queing jobs
                ctrl.cmd = 'task';
            end
            
            % =====================================================================
            % Load data For each frame
            % =====================================================================
            out_frm = {};
            retry_fields = {};
            
            for frm_id = 1:length(param.cmd.frms)
                
                % Create the array_proc output directories
                if ~exist(analysis_out_path,'dir')
                    mkdir(analysis_out_path);
                end
                
                % =====================================================================
                % Prepare task inputs
                % =====================================================================
                cur_frm = param.cmd.frms(frm_id);
                task_param.load.frm = cur_frm;
                
                file1=fullfile(ct_filename_out(param,param.analysis.in_path,''),sprintf('Data_%s_%03d.mat',param.day_seg,task_param.load.frm));
               
                schedulernotrun=0;   
                task_param = param;
                task_param.load.frm = cur_frm;
                if schedulernotrun~=0
                    % =====================================================================
                    % Setup the scheduler
                    % =====================================================================
                    
                    if strcmpi(param.sched.type,'custom_torque')
                        global ctrl; % Make this global for convenience in debugging
                        ctrl = torque_new_batch(param);
                        fprintf('Torque batch: %s\n', ctrl.batch_dir);
                        torque_compile('IceBedCoherenceIndex_tracker_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
                        
                    elseif ~strcmpi(param.sched.type,'no scheduler')
                        fd = [get_filenames(param.path,'','','.m',struct('recursive',1)); ...
                            get_filenames(param.path,'','','.mexa64',struct('recursive',1))];
                        fd_override = [get_filenames(param.path_override,'','','.m',struct('recursive',1)); ...
                            get_filenames(param.path_override,'','','.mexa64',struct('recursive',1))];
                        
                        fd = merge_filelists(fd, fd_override);
                        
                        % Remove SVN files from path
                        non_svn_path = zeros(size(fd));
                        for fn_idx = 1:length(fd)
                            if isempty(strfind(fd{fn_idx},'.svn'))
                                non_svn_path(fn_idx) = 1;
                            end
                        end
                        fd = fd(non_svn_path);
                        
                        % Initialize submission ctrl structure
                        global ctrl; % Make this global for convenience in debugging
                        ctrl = [];
                        ctrl.cmd = 'init';
                        ctrl.sched = param.sched;
                        ctrl.fd = fd;
                        ctrl = create_task(ctrl);
                        
                        % Prepare submission ctrl structure for queing jobs
                        ctrl.cmd = 'task';
                    end
                    
                    out_frm = {};
                    retry_fields = {};
                    
                    
                    % Create the array_proc output directories
                    if ~exist(analysis_out_path,'dir')
                        mkdir(analysis_out_path);
                    end
                    
                    % =====================================================================
                    % Prepare task inputs
                    % =====================================================================
                    cur_frm = param.cmd.frms(frm_id);
                    task_param.load.frm = cur_frm;
                end
                
                
                % =================================================================
                % Rerun only mode: Test to see if we need to run this task
                if param.sched.rerun_only
                    % If we are in rerun only mode AND all the IceBedCoherenceIndex_tracker task output files
                    % already exists, then we do not run the task
                    file_exists = true;
                    
                    if param.analysis.IceBedCoherenceIndex.en
                        out_fn = fullfile(ct_filename_out(param, ...
                            param.analysis.out_path, 'CSARP_basal_condition'), ...
                            sprintf('IceBedCoherenceIndex_%s_%03d.mat',param.day_seg,cur_frm));
                    end
                    
                    if ~exist(out_fn,'file')
                        file_exists = false;
                    end
                    
                    if file_exists
                        fprintf('IceBedCoherenceIndex_%s_%03d.mat- Already exists [rerun_only skipping] (%s)\n ',...
                            param.day_seg,cur_frm, datestr(now));
                        continue;
                    end
                end
                
                % =================================================================
                % Execute tasks/jobs
                
                fh = @IceBedCoherenceIndex_tracker_task;
                arg{1} = task_param;
                
                if strcmp(param.sched.type,'custom_torque')
                    create_task_param.conforming = true;
                    create_task_param.notes = sprintf('%s frames %d of %d total', ...
                        param.day_seg, cur_frm, length(frames.frame_idxs));
                    ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
                    
                elseif ~strcmp(param.sched.type,'no scheduler')
                    [ctrl,job_id,task_id] = create_task(ctrl,fh,1,arg);
                    fprintf('%s frames %d of %d total in job,task %d,%d (%s)\n', ...
                        param.day_seg, cur_frm, length(param.cmd.frms), job_id, task_id, datestr(now));
                    retry_fields{job_id,task_id}.arg = arg;
                    out_frm{end + 1} = cur_frm;
                    retry_fields{job_id,task_id}.out_idx = length(out_frm);
                else
                    fprintf('\n  %s frames  %d of %d total(%s)\n', ...
                        param.day_seg, cur_frm, length(param.cmd.frms), datestr(now));
                    [success] = fh(arg{1});
                end
            end
            
            % =======================================================================
            % Wait for jobs to complete if a scheduler was used
            % =======================================================================
            if strcmpi(param.sched.type,'custom_torque')
                % Wait until all submitted jobs to complete
                ctrl = torque_rerun(ctrl);
                if ~all(ctrl.error_mask == 0)
                    if ctrl.sched.stop_on_fail
                        torque_cleanup(ctrl);
                        error('Not all jobs completed, but out of retries (%s)', datestr(now));
                    else
                        warning('Not all jobs completed, but out of retries (%s)', datestr(now));
                        keyboard;
                    end
                else
                    fprintf('Jobs completed (%s)\n\n', datestr(now));
                end
                torque_cleanup(ctrl);
                
            elseif ~strcmpi(param.sched.type,'no scheduler')
                % ======================================================================
                % Wait for jobs to finish and clean up
                ctrl.cmd = 'done';
                ctrl = create_task(ctrl);
                if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
                    % Quit if a bad error occurred
                    fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
                    if strcmp(ctrl.sched.type,'torque')
                        fprintf('Often on the Torque scheduler, these are not bad errors\n');
                        fprintf('because of system instabilities (e.g. file IO failure)\n');
                        fprintf('and the task simply needs to be resubmitted. If this is the case,\n');
                        fprintf('run "ctrl.error_mask = 2" and then run "dbcont".\n');
                        keyboard
                        if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
                            return;
                        end
                    else
                        return
                    end
                end
                
                % ======================================================================
                % Retry jobs that failed
                retry = 1;
                while ctrl.error_mask == 2 && retry <= param.sched.max_retries
                    fprintf('Tasks failed, retry %d of max %d\n', retry, param.sched.max_retries);
                    
                    % Bookkeeping (move previous run information to "old_" variables)
                    old_ctrl = ctrl;
                    old_retry_fields = retry_fields;
                    retry_fields = {};
                    old_out_frms = out_frms;
                    out_frms = {};
                    
                    % Initialize submission ctrl structure
                    ctrl = [];
                    ctrl.cmd = 'init';
                    ctrl.sched = param.sched;
                    ctrl.fd = fd;
                    ctrl = create_task(ctrl);
                    
                    % Prepare submission ctrl structure for queing jobs
                    ctrl.cmd = 'task';
                    
                    % Submit failed tasks, but keep track of these in case they fail again
                    for job_idx = 1:length(old_ctrl.jobs)
                        for task_idx = old_ctrl.jobs{job_idx}.error_idxs
                            [ctrl,job_id,task_id] = create_task(ctrl,fh,2,old_retry_fields{job_idx,task_idx}.arg);
                            out_idx = old_retry_fields{job_idx,task_idx}.out_idx;
                            fprintf('  %d/%d: Processing records %d to %d in job,task %d,%d (%s)\n', ...
                                old_retry_fields{job_idx,task_idx}.frm, old_retry_fields{job_idx,task_idx}.break_idx, ...
                                old_out_frms{out_idx}(1), old_out_frms{out_idx}(end), ...
                                job_id, task_id, datestr(now));
                            retry_fields{job_id,task_id} = old_retry_fields{job_idx,task_idx};
                            out_frms{end + 1} = old_out_frms{out_idx};
                            retry_fields{job_id,task_id}.out_idx = length(out_frms);
                        end
                    end
                    
                    % Wait for tasks to complete and then cleanup
                    ctrl.cmd = 'done';
                    ctrl = create_task(ctrl);
                    retry = retry + 1;
                    
                    if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
                        % Quit if a bad error occurred
                        fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
                        if strcmp(ctrl.sched.type,'torque')
                            fprintf('Often on the Torque scheduler, these are not bad errors\n');
                            fprintf('because of system instabilities (e.g. file IO failure)\n');
                            fprintf('and the task simply needs to be resubmitted. If this is the case,\n');
                            fprintf('run "ctrl.error_mask = 2" and then run "dbcont".\n');
                            keyboard
                            if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
                                return;
                            end
                        else
                            return
                        end
                    end
                end
                if ctrl.error_mask ~= 0
                    fprintf('Not all jobs completed, but out of retries (%s)\n', datestr(now));
                    return;
                else
                    fprintf('Jobs completed (%s)\n\n', datestr(now));
                end
            end
        end
    end
    % concatnate coherence points from all frames of a segment
    
    if param.analysis.collate_en(1)
      if ~exist('bt_coherenceindex','var')
        bt_coherenceindex.GPS_time = [];
        bt_coherenceindex.Latitude = [];
        bt_coherenceindex.Longitude = [];
        bt_coherenceindex.Elevation_new = [];
        bt_coherenceindex.Time_new = [];
        bt_coherenceindex.Surface_new = [];
        bt_coherenceindex.SurfaceTime_new = [];
        bt_coherenceindex.Bottom_new = [];
        bt_coherenceindex.BottomTime_new = [];
        bt_coherenceindex.value = [];
        bt_coherenceindex.Abruptiveindex=[];
        bt_coherenceindex.AdjustedIntensity=[];
        bt_coherenceindex.Longitude_mean=[];
        bt_coherenceindex.Latitude_mean=[];
        bt_coherenceindex.frm_id = {};       
      end
      out_fn = fullfile(ct_filename_out(param,param.analysis.out_path, 'CSARP_basal_condition'), ...
        sprintf('IceBedCoherenceIndex_%s.mat',param.day_seg));
      out_fn = regexprep(out_fn,['/',param.day_seg],'');
       if exist(out_fn,'file')
           continue;
       end
      for frm_id = 1:length(param.cmd.frms)
        
        cur_frm = param.cmd.frms(frm_id);
        
        in_fn=  fullfile(ct_filename_out(param, ...
          param.analysis.out_path, 'CSARP_basal_condition'), ...
          sprintf('IceBedCoherenceIndex_%s_%03d.mat',param.day_seg,cur_frm));
        
        file_not_exists = false;
        
        if ~exist(in_fn,'file')
          file_not_exists = true;
        end
        
        if file_not_exists
          fprintf('IceBedCoherenceIndex_%s_%03d.mat- does not exists - run IceBedCoherenceIndex_analysis first (%s)\n ',...
            param.day_seg,cur_frm, datestr(now));
          continue;
        end
        load(in_fn) ;
        bt_coherenceindex.GPS_time = cat(2, bt_coherenceindex.GPS_time, bt.coherenceindex.GPS_time);
        bt_coherenceindex.Longitude = cat(2, bt_coherenceindex.Longitude, bt.coherenceindex.Longitude);
        bt_coherenceindex.Latitude = cat(2, bt_coherenceindex.Latitude, bt.coherenceindex.Latitude);
        bt_coherenceindex.Elevation_new = cat(2, bt_coherenceindex.Elevation_new, bt.coherenceindex.Elevation_new);
        bt_coherenceindex.Time_new = cat(2, bt_coherenceindex.Time_new, bt.coherenceindex.Time_new);
        bt_coherenceindex.Surface_new = cat(2, bt_coherenceindex.Surface_new, bt.coherenceindex.Surface_new);
        bt_coherenceindex.SurfaceTime_new = cat(2, bt_coherenceindex.SurfaceTime_new, bt.coherenceindex.SurfaceTime_new);
        bt_coherenceindex.Bottom_new = cat(2, bt_coherenceindex.Bottom_new, bt.coherenceindex.Bottom_new);
        bt_coherenceindex.BottomTime_new = cat(2, bt_coherenceindex.BottomTime_new, bt.coherenceindex.BottomTime_new);
        bt_coherenceindex.value = cat(2, bt_coherenceindex.value, bt.coherenceindex.value);
        bt_coherenceindex.Abruptiveindex = cat(2, bt_coherenceindex.Abruptiveindex, bt.coherenceindex.Abruptiveindex);
         bt_coherenceindex.AdjustedIntensity = cat(2, bt_coherenceindex.AdjustedIntensity, bt.coherenceindex.AdjustedIntensity);
         bt_coherenceindex.Latitude_mean = cat(2, bt_coherenceindex.Latitude_mean, bt.coherenceindex.Latitude_mean);
        bt_coherenceindex.Longitude_mean = cat(2, bt_coherenceindex.Longitude_mean, bt.coherenceindex.Longitude_mean);
        
    
        if ~isempty(bt_coherenceindex.value)
           for idx = 1:length(bt_coherenceindex.value)-length(bt_coherenceindex.frm_id)
                bt_coherenceindex.frm_id = cat(2,bt_coherenceindex.frm_id,sprintf('%s_%03d',param.day_seg,cur_frm));
           end
        end
      end
      bt.coherenceindex = bt_coherenceindex;
      param_analysis=params.analysis;
      fprintf('  Saving outputs %s\n', out_fn);  
      save(out_fn,'bt','param_analysis');  
      clear bt_coherence;
    end
  end
  
  % concatnate coherence points from all segments
  if param.analysis.collate_en(2)
    collate_en2_flag = 1;
    if ~exist('bt_all_coherenceindex','var')
        bt_all_coherenceindex.GPS_time = [];
        bt_all_coherenceindex.Latitude = [];
        bt_all_coherenceindex.Longitude = [];
        bt_all_coherenceindex.Elevation_new = [];
        bt_all_coherenceindex.Time_new = [];
        bt_all_coherenceindex.Surface_new = [];
        bt_all_coherenceindex.SurfaceTime_new = [];
        bt_all_coherenceindex.Bottom_new = [];
        bt_all_coherenceindex.BottomTime_new = [];
        bt_all_coherenceindex.value = [];
        bt_all_coherenceindex.frm_id = {};  
        bt_all_coherenceindex.GPS_time_ave = [];
     
      bt_all_coherenceindex.waveform = [];
       bt_all_coherenceindex.Latitude_mean=[];
        bt_all_coherenceindex.Longitude_mean=[];
         bt_all_coherenceindex.Abruptiveindex =[];
            bt_all_coherenceindex.AdjustedIntensity =[];
    end
    in_fn=  fullfile(ct_filename_out(param, param.analysis.out_path, 'CSARP_basal_condition'), ...
      sprintf('IceBedCoherenceIndex_%s.mat',param.day_seg));
    in_fn = regexprep(in_fn,['/',param.day_seg],'');
    
    file_not_exists = false;
    if ~exist(in_fn,'file')
      file_not_exists = true;
    end
    
    if file_not_exists
      fprintf('IceBedCoherenceIndex_%s.mat- does not exists - run IceBedCoherenceIndex_analysis first (%s)\n ',...
        param.day_seg, datestr(now));
      continue;
    end
    load(in_fn) ;
        bt_all_coherenceindex.GPS_time = cat(2, bt_all_coherenceindex.GPS_time, bt.coherenceindex.GPS_time);
         
        bt_all_coherenceindex.Longitude = cat(2, bt_all_coherenceindex.Longitude, bt.coherenceindex.Longitude);
        bt_all_coherenceindex.Latitude = cat(2, bt_all_coherenceindex.Latitude, bt.coherenceindex.Latitude);
        bt_all_coherenceindex.Elevation_new = cat(2, bt_all_coherenceindex.Elevation_new, bt.coherenceindex.Elevation_new);
        bt_all_coherenceindex.Time_new = cat(2, bt_all_coherenceindex.Time_new, bt.coherenceindex.Time_new);
        bt_all_coherenceindex.Surface_new = cat(2, bt_all_coherenceindex.Surface_new, bt.coherenceindex.Surface_new);
        bt_all_coherenceindex.SurfaceTime_new = cat(2, bt_all_coherenceindex.SurfaceTime_new, bt.coherenceindex.SurfaceTime_new);
        bt_all_coherenceindex.Bottom_new = cat(2, bt_all_coherenceindex.Bottom_new, bt.coherenceindex.Bottom_new);
        bt_all_coherenceindex.BottomTime_new = cat(2, bt_all_coherenceindex.BottomTime_new, bt.coherenceindex.BottomTime_new);
       bt_all_coherenceindex.value = cat(2, bt_all_coherenceindex.value, bt.coherenceindex.value);
         bt_all_coherenceindex.Latitude_mean = cat(2, bt_all_coherenceindex.Latitude_mean, bt.coherenceindex.Latitude_mean);
        bt_all_coherenceindex.Longitude_mean = cat(2, bt_all_coherenceindex.Longitude_mean, bt.coherenceindex.Longitude_mean);
           bt_all_coherenceindex.Abruptiveindex = cat(2, bt_all_coherenceindex.Abruptiveindex, bt.coherenceindex.Abruptiveindex);
           bt_all_coherenceindex.AdjustedIntensity = cat(2, bt_all_coherenceindex.AdjustedIntensity, bt.coherenceindex.AdjustedIntensity);
  
  end
end

if exist('collate_en2_flag','var')
  if collate_en2_flag
    out_fn = fullfile(ct_filename_out(param,param.analysis.out_path, 'CSARP_basal_condition'), ...
      sprintf('IceBedCoherenceIndex_%s.mat',param.season_name));
    out_fn = regexprep(out_fn,['/',param.day_seg],'');
    clear bt;
    bt.coherenceindex = bt_all_coherenceindex;
    fprintf('  Saving outputs %s\n', out_fn);
    param_analysis = params.analysis;
    save(out_fn,'bt','param_analysis');
    clear bt_all_coherenceindex;
  end
end

if param.analysis.generate_map.en 
  geotiff_fn = ct_filename_gis(param, param.analysis.map_path);
  proj = geotiffinfo(geotiff_fn);
  [A CMAP R]= geotiffread(geotiff_fn);
  figure
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  hold on;
  [gps.x,gps.y] = projfwd(proj,bt.coherenceindex.Latitude_mean,bt.coherenceindex.Longitude_mean);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  grid
  title('Coherenceindex');
  xlabel('Latitude(deg)');
  ylabel('Longitude(deg)');
  hold on;
  scatter(gps.x,gps.y,20,bt.coherenceindex.value,'fill')
  colorbar
  if length(bt.coherenceindex.value) >1
    caxis([min(bt.coherenceindex.value),max(bt.coherenceindex.value)]);
  elseif length(bt.coherenceindex.value) == 1
    caxis([bt.coherenceindex.value-50,bt.coherenceindex.value+50]);
  end
  
  %abr
  figure
  mapshow(rgb2gray(A),CMAP/1e3);
  xlabel('X (km)');
  ylabel('Y (km)');
  hold on;
  [gps.x,gps.y] = projfwd(proj,bt.coherenceindex.Latitude_mean,bt.coherenceindex.Longitude_mean);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  grid
  title('Abruptive index');
  xlabel('Latitude(deg)');
  ylabel('Longitude(deg)');
  hold on;
  scatter(gps.x,gps.y,20,bt.coherenceindex.Abruptiveindex,'fill')
  colorbar
  if length(bt.coherenceindex.Abruptiveindex) >1
    caxis([min(bt.coherenceindex.Abruptiveindex),max(bt.coherenceindex.Abruptiveindex)]);
  elseif length(bt.coherenceindex.Abruptiveindex) == 1
    caxis([bt.coherenceindex.Abruptiveindex-50,bt.coherenceindex.Abruptiveindex+50]);
  end
end

return;



