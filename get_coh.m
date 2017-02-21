

function get_coh(params,param_override)

global gRadar;

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

        
        % =====================================================================
        % Setup processing
        % =====================================================================
        
        % Get WGS84 ellipsoid parameters
        physical_constants;
        
        if ~isfield(param.records,'records_fn')
            param.records.records_fn = '';
        end
        if ~isfield(param.records,'frames_fn')
            param.records.frames_fn = '';
        end
        
        % Load frames file
        load(ct_filename_support(param,param.records.frames_fn,'frames'));
        % Load records file
        records_fn = ct_filename_support(param,param.records.records_fn,'records');
        records = load(records_fn);
        
        global g_data;
        g_data = [];
        
        qlook_out_path = ct_filename_out(param, param.analysis.out_path);
        
        if isempty(param.cmd.frms)
            param.cmd.frms = 1:length(frames.frame_idxs);
        end
        % Remove frames that do not exist from param.cmd.frms list
        [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
        if length(valid_frms) ~= length(param.cmd.frms)
            bad_mask = ones(size(param.cmd.frms));
            bad_mask(keep_idxs) = 0;
            warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
                param.cmd.frms(find(bad_mask,1)));
            param.cmd.frms = valid_frms;
        end
        
        if ~isfield(param.sched,'rerun_only') || isempty(param.sched.rerun_only)
            param.sched.rerun_only = false;
        end
        
        if ~isfield(param.get_heights,'ground_based')
            param.get_heights.ground_based = [];
        end
        
        if ~isfield(param.get_heights.qlook,'save_format') || isempty(param.get_heights.qlook.save_format)
            param.get_heights.qlook.save_format = '7.3';
        end
        save_format = sprintf('-v%s',param.get_heights.qlook.save_format);
        
        if isfield(param.get_heights,'deconvolution') ...
                && ~isempty(param.get_heights.deconvolution) ...
                && param.get_heights.deconvolution == 3
            %% Get version information out of the deconvolution file
            out_fn_dir = ct_filename_out(param,'', 'CSARP_noise');
            out_segment_fn_dir = fileparts(out_fn_dir);
            out_segment_fn = fullfile(out_segment_fn_dir,sprintf('deconv_%s.mat', param.day_seg));
            spec = load(out_segment_fn,'param_collate');
            
            param.get_heights.deconvolution_sw_version = spec.param_collate.sw_version;
            param.get_heights.deconvolution_params = spec.param_collate.analysis.specular;
        end
        
        if isfield(param.get_heights,'coh_noise_method') ...
                && ~isempty(param.get_heights.coh_noise_method) ...
                && any(param.get_heights.coh_noise_method == [17 19])
            %% Get version information out of the coherent noise file
            
            cdf_fn_dir = fileparts(ct_filename_out(param,param.get_heights.coh_noise_arg{4}, ''));
            cdf_fn = fullfile(cdf_fn_dir,sprintf('coh_noise_simp_%s.nc', param.day_seg));
            
            tmp = netcdf_to_mat(cdf_fn,[],'^sw_version.*');
            param.get_heights.coh_noise_version = tmp.sw_version;
            tmp = netcdf_to_mat(cdf_fn,[],'^param_collate.*');
            param.get_heights.coh_noise_params = tmp.param_collate;
        end
        
        % Cleanup folders
        if ~param.sched.rerun_only
            if exist(qlook_out_path,'dir')
                for frm = param.cmd.frms
                    del_paths = get_filenames(qlook_out_path,sprintf('ql_data_%03d',frm),'','',struct('type','d'));
                    for idx = 1:length(del_paths)
                        fprintf('Removing path: %s\n', del_paths{idx});
                        rmdir(del_paths{idx},'s');
                    end
                end
            end
        end
        
        seg_file = fullfile(ct_filename_out(param,param.analysis.out_path, 'CSARP_basal_condition'), ...
            sprintf('IceBedCoherenceIndex_%s.mat',param.day_seg));
        seg_file= regexprep(seg_file,['/',param.day_seg],'');
        existdir=0;
        if exist(seg_file,'file')
            existdir = 1;
        end
        
        if existdir==0
            
            
            % =====================================================================
            % Setup static inputs for get_heights_task
            % =====================================================================
            task_param = param;
            task_param.load.imgs = param.get_heights.imgs;
            
            
            % =====================================================================
            % Setup the scheduler
            % =====================================================================
            
            if strcmpi(param.sched.type,'custom_torque')
                global ctrl; % Make this global for convenience in debugging
                ctrl = torque_new_batch(param);
                fprintf('Torque batch: %s\n', ctrl.batch_dir);
                torque_compile('get_coh_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
                
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
            
            if ~strcmpi(param.sched.type,'no scheduler')
                param.get_heights.surf.manual = 0; % Turn manual pick off
            end
            
            % =====================================================================
            % Load data and run get_heights tasks
            %
            % For each frame load REC_BLOCK_SIZE records at a time (code groups
            % by file index, but has to watch negative offset values which imply
            % the record starts in a previous file and carries over into the next)
            %    --> The last block can be up to 2*REC_BLOCK_SIZE
            % =====================================================================
            out_recs = {};
            retry_fields = {};
            brk_length=zeros(length(param.cmd.frms));
            
            for frm_idx = 1:length(param.cmd.frms)
                frm = param.cmd.frms(frm_idx);
                
                file1 = fullfile(ct_filename_out(param, ...
                    param.analysis.out_path,'CSARP_basal_condition'), ...
                    sprintf('IceBedCoherenceIndex_%s_%03d.mat',param.day_seg,frm));
                
                if exist(file1,'file')
                    success = true;
                    continue;
                end
               param.get_heights.frm_types= {0,[0 1],0,0,0};

                % Check digits of proc_mode from frames file and make sure the user has
                % specified to process this frame type
                if ct_proc_frame(frames.proc_mode(frm),param.get_heights.frm_types)
                    fprintf('get_heights %s_%03i (%i of %i) %s\n', param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
                else
                    fprintf('Skipping %s_%03i (no process frame)\n', param.day_seg, frm);
                    continue;
                end
                
                
                
                if frm < length(frames.frame_idxs)
                    recs = frames.frame_idxs(frm):frames.frame_idxs(frm+1);
                else
                    recs = frames.frame_idxs(frm):length(records.lat);
                end
                
                % Determine where breaks in processing are going to occur
                rec = recs(1);
                
                if isempty(param.get_heights.block_size)
                    REC_BLOCK_SIZE = 10000;
                else
                    if numel(param.get_heights.block_size) == 2
                        REC_BLOCK_SIZE = param.get_heights.block_size(1);
                        block_overlap = param.get_heights.block_size(2);
                    else
                        REC_BLOCK_SIZE = param.get_heights.block_size;
                        block_overlap = 0;
                    end
                end
                if length(recs) < 2*REC_BLOCK_SIZE
                    breaks = 1;
                else
                    breaks = 1:REC_BLOCK_SIZE:length(recs)-REC_BLOCK_SIZE;
                end
                
                task_param.proc.frm = frm;
                
                % Begin loading data
                for break_idx = 1:length(breaks)
                    % Determine the current records being processed
                    if break_idx < length(breaks)
                        cur_recs_keep = [recs(breaks(break_idx)) recs(breaks(break_idx+1)-1)];
                        cur_recs = [max(1,recs(breaks(break_idx))-block_overlap) ...
                            recs(breaks(break_idx+1)-1)+block_overlap];
                    else
                        cur_recs_keep = [recs(breaks(break_idx)) recs(end)];
                        cur_recs = [max(1,recs(breaks(break_idx))-block_overlap) recs(end)];
                    end
                    
                    % =====================================================================
                    % Prepare task inputs
                    % =====================================================================
                    task_param.load.recs = cur_recs;
                    task_param.load.recs_keep = cur_recs_keep;
                    task_param.load.blockno=break_idx;
                    task_param.load.imgs={[2 9; 2 10; 2 11; 2 12]};
                    task_param.get_heights.imgs={[2 9; 2 10; 2 11; 2 12]};
                    task_param.get_heights.B_filter=hanning(11).'/sum(hanning(11));
                    task_param.get_heights.decimate_factor=7;
                    task_param.get_heights.inc_ave=0;
                    task_param.get_heights.sur.en=0;
                    
                    
                    file1 = fullfile(ct_filename_out(param, ...
                        param.analysis.out_path,'CSARP_basal_condition'), ...
                        sprintf('IceBedCoherenceIndex_%s_%03d_%02d.mat',param.day_seg,task_param.proc.frm,task_param.load.blockno));
                    
                    
                    if exist(file1,'file')
                        success = true;
                        continue;
                    end
                    
                    
                    % =================================================================
                    % Execute tasks/jobs
                    fh = @get_coh_task;
                    if isfield(frames,'nyquist_zone') && ~isnan(frames.nyquist_zone(frm))
                        task_param.radar.wfs(1).nyquist_zone = frames.nyquist_zone(frm);
                    elseif isfield(param.radar.wfs(1),'nyquist_zone') ...
                            && ~isempty(param.radar.wfs(1).nyquist_zone) ...
                            && ~isnan(param.radar.wfs(1).nyquist_zone)
                        task_param.radar.wfs(1).nyquist_zone = param.radar.wfs(1).nyquist_zone;
                    end
                    arg{1} = task_param;
                    
                    if strcmp(param.sched.type,'custom_torque')
                        create_task_param.conforming = true;
                        create_task_param.notes = sprintf('%s_%03d (%d of %d)/%d of %d records %d-%d', ...
                            param.day_seg, frm, frm_idx, length(param.cmd.frms), break_idx, length(breaks), cur_recs(1), cur_recs(end));
                        ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
                        
                    elseif ~strcmp(param.sched.type,'no scheduler')
                        [ctrl,job_id,task_id] = create_task(ctrl,fh,1,arg);
                        fprintf('  %d/%d: records %d to %d in job,task %d,%d (%s)\n', ...
                            frm, break_idx, cur_recs(1), cur_recs(end), job_id, task_id, datestr(now));
                        retry_fields{job_id,task_id}.frm = frm;
                        retry_fields{job_id,task_id}.break_idx = break_idx;
                        retry_fields{job_id,task_id}.arg = arg;
                        out_recs{end + 1} = cur_recs;
                        retry_fields{job_id,task_id}.out_idx = length(out_recs);
                    else
                        fprintf('  %s_%03d (%d of %d)/%d of %d: records %d-%d (%s)\n', ...
                            param.day_seg, frm, frm_idx, length(param.cmd.frms), break_idx, length(breaks), cur_recs(1), cur_recs(end), datestr(now));
                        [success] = fh(arg{1});
                    end
                    
                end
                brk_length(frm_idx)= length(breaks);
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
                    old_out_recs = out_recs;
                    out_recs = {};
                    
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
                                old_out_recs{out_idx}(1), old_out_recs{out_idx}(end), ...
                                job_id, task_id, datestr(now));
                            retry_fields{job_id,task_id} = old_retry_fields{job_idx,task_idx};
                            out_recs{end + 1} = old_out_recs{out_idx};
                            retry_fields{job_id,task_id}.out_idx = length(out_recs);
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

            
            %%  %Concatenate all blocks into a frame
            for frm_idx = 1:length(param.cmd.frms)
                for i=1:brk_length(frm_idx)
                    if ~exist('bt_coherenceindex','var')
                        bt_coherenceindex.GPS_time = [];
                        bt_coherenceindex.Latitude = [];
                        bt_coherenceindex.Longitude = [];
                        bt_coherenceindex.Elevation_new = [];
                        
                        bt_coherenceindex.Surface_new = [];
                        bt_coherenceindex.SurfaceTime_new = [];
                        bt_coherenceindex.Bottom_new = [];
                        bt_coherenceindex.BottomTime_new = [];
                        bt_coherenceindex.value = [];
                        bt_coherenceindex.Abruptiveindex=[];
                        bt_coherenceindex.AdjustedIntensity=[];
                        bt_coherenceindex.Longitude_mean=[];
                        bt_coherenceindex.Latitude_mean=[];
                        
                    end
                    if ~exist('bottom','var')
                        bottom.peakval = [];
                        bottom.peakidx = [];
                        bottom.waveform = [];
                        bottom.inc_wf_ave= [];
                        bottom.avepeakindex= [];
                    end
                    
                    out_fn = fullfile(ct_filename_out(param,param.analysis.out_path, 'CSARP_basal_condition'), ...
                        sprintf('IceBedCoherenceIndex_%s_%03d.mat',param.day_seg, param.cmd.frms(frm_idx)));
                    out_fn = regexprep(out_fn,['/',param.day_seg,param.cmd.frms(frm_idx)],'');
                    
                    
                    
                    in_fn=  fullfile(ct_filename_out(param, ...
                        param.analysis.out_path, 'CSARP_basal_condition'), ...
                        sprintf('IceBedCoherenceIndex_%s_%03d_%02d.mat',param.day_seg, param.cmd.frms(frm_idx),i));
                    
                    file_not_exists = false;
                    
                    if ~exist(in_fn,'file')
                        file_not_exists = true;
                    end
                    
                    if file_not_exists
                        fprintf('IceBedCoherenceIndex_%s_%03d%02d.mat- does not exists - run IceBedCoherenceIndex_analysis first (%s)\n ',...
                            param.day_seg,param.cmd.frms(frm_idx), datestr(now));
                        continue;
                    end
                    load(in_fn) ;
                    
                    
                    bt_coherenceindex.GPS_time = cat(2, bt_coherenceindex.GPS_time, bt.coherenceindex.GPS_time);
                    bt_coherenceindex.Longitude = cat(2, bt_coherenceindex.Longitude, bt.coherenceindex.Longitude);
                    bt_coherenceindex.Latitude = cat(2, bt_coherenceindex.Latitude, bt.coherenceindex.Latitude);
                    bt_coherenceindex.Elevation_new = cat(2, bt_coherenceindex.Elevation_new, bt.coherenceindex.Elevation_new);
                    bt_coherenceindex.Time_new = bt.coherenceindex.Time_new;
                    bt_coherenceindex.Surface_new = cat(2, bt_coherenceindex.Surface_new, bt.coherenceindex.Surface_new);
                    bt_coherenceindex.SurfaceTime_new = cat(2, bt_coherenceindex.SurfaceTime_new, bt.coherenceindex.SurfaceTime_new);
                    bt_coherenceindex.Bottom_new = cat(2, bt_coherenceindex.Bottom_new, bt.coherenceindex.Bottom_new);
                    bt_coherenceindex.BottomTime_new = cat(2, bt_coherenceindex.BottomTime_new, bt.coherenceindex.BottomTime_new);
                    bt_coherenceindex.value = cat(2, bt_coherenceindex.value, bt.coherenceindex.value);
                    bt_coherenceindex.Abruptiveindex = cat(2, bt_coherenceindex.Abruptiveindex, bt.coherenceindex.Abruptiveindex);
                     bt_coherenceindex.AdjustedIntensity = cat(2, bt_coherenceindex.AdjustedIntensity, bt.coherenceindex.AdjustedIntensity);
                    bt_coherenceindex.Latitude_mean = cat(2, bt_coherenceindex.Latitude_mean, bt.coherenceindex.Latitude_mean);
                    bt_coherenceindex.Longitude_mean = cat(2, bt_coherenceindex.Longitude_mean, bt.coherenceindex.Longitude_mean);
                    
                    bottom.peakval = cat(2,bottom.peakval,bt.peakval);
                    bottom.peakidx = cat(2,bottom.peakidx,bt.peakidx);
                    bottom.waveform = cat(2,bottom.waveform,bt.waveform);
                    bottom.inc_wf_ave= cat(2,bottom.inc_wf_ave,bt.inc_wf_ave);
                    bottom.avepeakindex= cat(2,bottom.avepeakindex,bt.avepeakindex);
                    
                    
                    bt=bottom;
                    bt.coherenceindex = bt_coherenceindex;
                    param_analysis=params.analysis;
                    
                    if i==brk_length(frm_idx)
                        out_fn = fullfile(ct_filename_out(param,param.analysis.out_path, 'CSARP_basal_condition'), ...
                            sprintf('IceBedCoherenceIndex_%s_%03d.mat',param.day_seg, param.cmd.frms(frm_idx)));
                        
                        
                        fprintf('  Saving outputs %s\n', out_fn);
                        save(out_fn,'bt','param_analysis');
                        clear bt_coherenceindex;
                        clear bottom;
                        for j=1:brk_length(frm_idx)
                            in_fn=  fullfile(ct_filename_out(param, ...
                                param.analysis.out_path, 'CSARP_basal_condition'), ...
                                sprintf('IceBedCoherenceIndex_%s_%03d_%02d.mat',param.day_seg, param.cmd.frms(frm_idx),j));
                            delete(in_fn);
                        end
                    end
                end
            end
     end
 
     %% concatnate coherence points from all frames of a segment
     
     if param.analysis.collate_en(1)
         if ~exist('bt_coherenceindex','var')
             %  bt_coherenceindex.GPS_time = [];
             % bt_coherenceindex.Latitude = [];
             % bt_coherenceindex.Longitude = [];
             % bt_coherenceindex.Elevation_new = [];
             
             %  bt_coherenceindex.Surface_new = [];
             %  bt_coherenceindex.SurfaceTime_new = [];
             % bt_coherenceindex.Bottom_new = [];
             % bt_coherenceindex.BottomTime_new = [];
             bt_coherenceindex.value = [];
             bt_coherenceindex.Abruptiveindex=[];
             bt_coherenceindex.AdjustedIntensity=[];
             bt_coherenceindex.Longitude_mean=[];
             bt_coherenceindex.Latitude_mean=[];
             
         end
         %       if ~exist('bottom','var')
         %              bottom.peakval = [];
         %           bottom.peakidx = [];
         %            bottom.waveform = [];
         %             bottom.inc_wf_ave= [];
         %              bottom.avepeakindex= [];
         %       end
         
         out_fn = fullfile(ct_filename_out(param,param.analysis.out_path, 'CSARP_basal_condition'), ...
             sprintf('IceBedCoherenceIndex_%s.mat',param.day_seg));
         out_fn = regexprep(out_fn,['/',param.day_seg],'');
         for frm_id = 1:length(frames.frame_idxs)
             
             cur_frm =frm_id;
             
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
             %          bt_coherenceindex.GPS_time = cat(2, bt_coherenceindex.GPS_time, bt.coherenceindex.GPS_time);
             %         bt_coherenceindex.Longitude = cat(2, bt_coherenceindex.Longitude, bt.coherenceindex.Longitude);
             %         bt_coherenceindex.Latitude = cat(2, bt_coherenceindex.Latitude, bt.coherenceindex.Latitude);
             %         bt_coherenceindex.Elevation_new = cat(2, bt_coherenceindex.Elevation_new, bt.coherenceindex.Elevation_new);
             %         bt_coherenceindex.Time_new = bt.coherenceindex.Time_new;
             %         bt_coherenceindex.Surface_new = cat(2, bt_coherenceindex.Surface_new, bt.coherenceindex.Surface_new);
             %         bt_coherenceindex.SurfaceTime_new = cat(2, bt_coherenceindex.SurfaceTime_new, bt.coherenceindex.SurfaceTime_new);
             %         bt_coherenceindex.Bottom_new = cat(2, bt_coherenceindex.Bottom_new, bt.coherenceindex.Bottom_new);
             %         bt_coherenceindex.BottomTime_new = cat(2, bt_coherenceindex.BottomTime_new, bt.coherenceindex.BottomTime_new);
             bt_coherenceindex.value = cat(2, bt_coherenceindex.value, bt.coherenceindex.value);
             bt_coherenceindex.Abruptiveindex = cat(2, bt_coherenceindex.Abruptiveindex, bt.coherenceindex.Abruptiveindex);
             bt_coherenceindex.AdjustedIntensity = cat(2, bt_coherenceindex.AdjustedIntensity, bt.coherenceindex.AdjustedIntensity);
             bt_coherenceindex.Latitude_mean = cat(2, bt_coherenceindex.Latitude_mean, bt.coherenceindex.Latitude_mean);
             bt_coherenceindex.Longitude_mean = cat(2, bt_coherenceindex.Longitude_mean, bt.coherenceindex.Longitude_mean);
             %
             %          bottom.peakval = cat(2,bottom.peakval,bt.peakval);
             %           bottom.peakidx = cat(2,bottom.peakidx,bt.peakidx);
             %            bottom.waveform = cat(2,bottom.waveform,bt.waveform);
             %             bottom.inc_wf_ave= cat(2,bottom.inc_wf_ave,bt.inc_wf_ave);
             %              bottom.avepeakindex= cat(2,bottom.avepeakindex,bt.avepeakindex);
             
         end
         clear bt;
         %       bt=bottom;Padj(Padj<(max(Padj)-40))=max(Padj)-40;
   %   bt_coherenceindex.AdjustedIntensity( bt_coherenceindex.AdjustedIntensity<(max( bt_coherenceindex.AdjustedIntensity)-40))=max( bt_coherenceindex.AdjustedIntensity)-40;
         bt.coherenceindex = bt_coherenceindex;
         param_analysis=params.analysis;
         fprintf('  Saving outputs %s\n', out_fn);
         save(out_fn,'bt','param_analysis');
         clear bt_coherenceindex;
         clear bottom;
     end
    end
 end
 %% % concatnate coherence points from all segments
 if param.analysis.collate_en(2)
     collate_en2_flag = 1;
     
     out_fn = fullfile(ct_filename_out(param,param.analysis.out_path, 'CSARP_basal_condition'), ...
         sprintf('IceBedCoherenceIndex_%s.mat',param.season_name));
     out_fn = regexprep(out_fn,['/',param.day_seg],'');
     
     file_not_exists = false;
     if ~exist(out_fn,'file')
         file_not_exists = true;
     end
     
     if ~file_not_exists
         fprintf('IceBedCoherenceIndex_%s.mat- does not exists - run IceBedCoherenceIndex_analysis first (%s)\n ',...
             param.day_seg, datestr(now));
         continue;
     end
     
    
     
     if ~exist('bt_all_coherenceindex','var')
         %         bt_all_coherenceindex.GPS_time = [];
         %         bt_all_coherenceindex.Latitude = [];
         %         bt_all_coherenceindex.Longitude = [];
         %         bt_all_coherenceindex.Elevation_new = [];
         %         bt_all_coherenceindex.Surface_new = [];
         %         bt_all_coherenceindex.SurfaceTime_new = [];
         %         bt_all_coherenceindex.Bottom_new = [];
         %         bt_all_coherenceindex.BottomTime_new = [];
         bt_all_coherenceindex.value = [];
         bt_all_coherenceindex.Latitude_mean=[];
         bt_all_coherenceindex.Longitude_mean=[];
         bt_all_coherenceindex.Abruptiveindex =[];
         bt_all_coherenceindex.AdjustedIntensity =[];
     end
     %      if ~exist('bottom_all','var')
     %              bottom_all.peakval = [];
     %           bottom_all.peakidx = [];
     %            bottom_all.waveform = [];
     %             bottom_all.inc_wf_ave= [];
     %              bottom_all.avepeakindex= [];
     %      end
     
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
     %
     %         bt_all_coherenceindex.GPS_time = cat(2, bt_all_coherenceindex.GPS_time, bt.coherenceindex.GPS_time);
     %         bt_all_coherenceindex.Longitude = cat(2, bt_all_coherenceindex.Longitude, bt.coherenceindex.Longitude);
     %         bt_all_coherenceindex.Latitude = cat(2, bt_all_coherenceindex.Latitude, bt.coherenceindex.Latitude);
     %         bt_all_coherenceindex.Elevation_new = cat(2, bt_all_coherenceindex.Elevation_new, bt.coherenceindex.Elevation_new);
     %         bt_all_coherenceindex.Time_new = bt.coherenceindex.Time_new;
     %         bt_all_coherenceindex.Surface_new = cat(2, bt_all_coherenceindex.Surface_new, bt.coherenceindex.Surface_new);
     %         bt_all_coherenceindex.SurfaceTime_new = cat(2, bt_all_coherenceindex.SurfaceTime_new, bt.coherenceindex.SurfaceTime_new);
     %         bt_all_coherenceindex.Bottom_new = cat(2, bt_all_coherenceindex.Bottom_new, bt.coherenceindex.Bottom_new);
     %         bt_all_coherenceindex.BottomTime_new = cat(2, bt_all_coherenceindex.BottomTime_new, bt.coherenceindex.BottomTime_new);
     bt_all_coherenceindex.value = cat(2, bt_all_coherenceindex.value, bt.coherenceindex.value);
     bt_all_coherenceindex.Abruptiveindex = cat(2, bt_all_coherenceindex.Abruptiveindex, bt.coherenceindex.Abruptiveindex);
     bt_all_coherenceindex.AdjustedIntensity = cat(2, bt_all_coherenceindex.AdjustedIntensity, bt.coherenceindex.AdjustedIntensity);
     bt_all_coherenceindex.Latitude_mean = cat(2, bt_all_coherenceindex.Latitude_mean, bt.coherenceindex.Latitude_mean);
     bt_all_coherenceindex.Longitude_mean = cat(2, bt_all_coherenceindex.Longitude_mean, bt.coherenceindex.Longitude_mean);
     
     %          bottom_all.peakval = cat(2,bottom_all.peakval,bt.peakval);
     %           bottom_all.peakidx = cat(2,bottom_all.peakidx,bt.peakidx);
     %            bottom_all.waveform = cat(2,bottom_all.waveform,bt.waveform);
     %             bottom_all.inc_wf_ave= cat(2,bottom_all.inc_wf_ave,bt.inc_wf_ave);
     %              bottom_all.avepeakindex= cat(2,bottom_all.avepeakindex,bt.avepeakindex);
     
 end
 
end
if exist('collate_en2_flag','var')
    if collate_en2_flag
        out_fn = fullfile(ct_filename_out(param,param.analysis.out_path, 'CSARP_basal_condition'), ...
            sprintf('IceBedCoherenceIndex_%s.mat',param.season_name));
        out_fn = regexprep(out_fn,['/',param.day_seg],'');
        clear bt;
        %     bt=bottom_all;
        bt.coherenceindex = bt_all_coherenceindex;
        fprintf('  Saving outputs %s\n', out_fn);
        param_analysis = params.analysis;
        save(out_fn,'bt','param_analysis');
        clear bottom_all;
    end
end

if param_idx==length(params)
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
end

return







