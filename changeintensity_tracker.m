
function changeintensity_tracker(params,param_override)


for param_idx = 1:length(params)
 
param = params(param_idx);
param = merge_structs(param, param_override);
    if isfield(param.cmd,'generic') && param.cmd.generic
    if param.analysis.IceBedCoherenceIndex.en
        param = params(param_idx);
        param = merge_structs(param, param_override);
    
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
        
     
        if ~isfield(param.get_heights.qlook,'save_format') || isempty(param.get_heights.qlook.save_format)
            param.get_heights.qlook.save_format = '7.3';
        end
        save_format = sprintf('-v%s',param.get_heights.qlook.save_format);
        
        
 
   
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
     bt.coherenceindex.AdjustedIntensity= bt.coherenceindex.AdjustedIntensity+param.analysis.adjint;
     save(in_fn,'bt','param_analysis')
    end
  end
end
return