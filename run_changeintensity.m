% script run_changeintensity
%
% Script for running changerunintensity_tracker: find coherence index
% of the Ice Bed to determine the smoothness of the Ice Bed
% 
%
% Authors: Jilu Li
%
% See also: master.m, IceBedCoherenceIndex_tracker.m IceBedCoherenceIndex_tracker_task.m,
tic
dbstop if error

% =====================================================================
% Debug Setup
% =====================================================================
% change the file_name of xls file and use generic field of command...
% ...spreadsheet to calculate the roughness parameters
params = read_param_xls(ct_filename_param('/users/manjish/Documents/MATLAB/rds_param_2012_Greenland_P333.xls'),'',{'analysis_IceBedCoherenceIndex','analysis'});

clear('param_override');
%param_override.sched.type='custom_torque';
 param_override.sched.type = 'no scheduler';
param_override.sched.rerun_only = true;
param_override.sched.cluster_size=inf;
param_override.sched.stop_on_fail=true;

param_override.sched.submit_arguments    = '-l nodes=1:ppn=1,walltime=15:00';
% Input checking
if ~exist('params','var')
  error('A struct array of parameters must be passed in\n');
end
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

changeintensity_tracker(params,param_override);

return
