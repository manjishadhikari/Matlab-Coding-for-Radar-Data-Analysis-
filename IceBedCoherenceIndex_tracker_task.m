function [success] = IceBedCoherenceIndex_tracker_task(param)
% [success] = IceBedCoherenceIndex_tracker_task(param)
%
% Cluster task for IceBedCoherenceIndex_tracker. Does the actual data loading
% and ice bed coherence index analysis to find the coherence index 
% to find out the smoothness of ice bed
%
% param = struct controlling the loading, processing, and coherenceindexity parameter generation
%
% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
%
% Author: Jilu Li
%
% See also IceBedCoherenceIndex_tracker.m

global g_data;
physical_constants;

%%
if ~isfield(param.get_heights.qlook,'save_format') || isempty(param.get_heights.qlook.save_format)
  param.get_heights.qlook.save_format = '7.3';
end
save_format = sprintf('-v%s',param.get_heights.qlook.save_format);

%% =====================================================================
% Determine the frames to load

if param.analysis.IceBedCoherenceIndex.en
    
  frm_fn = fullfile(ct_filename_out(param,param.analysis.in_path,''),sprintf('Data_%s_%03d.mat',param.day_seg,param.load.frm));
    
  if any([strfind(param.cmd.notes,'High-altitude'),strfind(param.cmd.notes,'Tomography')]) & ~exist(frm_fn,'file') 
    frm = load(fullfile(ct_filename_out(param,param.analysis.in_path,''),sprintf('Data_img_01_%s_%03d.mat',param.day_seg,param.load.frm)));
  else
    frm = load(frm_fn);    
  end
  layer = load(fullfile(ct_filename_out(param,'CSARP_post/layerData',''), ...
      sprintf('Data_%s_%03d.mat',param.day_seg,param.load.frm)));
  
distance=geodetic_to_along_track(frm.Latitude,frm.Longitude,frm.Elevation); %Horizontal Distance
sf = interp1(layer.GPS_time,layer.layerData{1}.value{2}.data,frm.GPS_time,'linear','extrap'); %Two way prop time to surface
bttime = interp1(layer.GPS_time, layer.layerData{2}.value{2}.data,frm.GPS_time,'linear','extrap');  %Two way prop. time to Bottom
bttime(bttime>frm.Time(end) | bttime<frm.Time(1)) = NaN; % remove pick errors exist
  
power=lp(frm.Data);
  
for rline=1:size(frm.Data,2)
   if isnan(bttime(rline))
       continue;
   end
    bt_idx_m = find(frm.Time>bttime(rline),1,'first'); 
    [~,bt_idx]=max(power(bt_idx_m-10:bt_idx_m+10,rline));
    dt=frm.Time(2)-frm.Time(1);
    bttime_n(rline)=(bttime(rline)-10*dt)+(bt_idx-1)*dt;
end
 debug_flag=0;
  if debug_flag==1
      figure(1); imagesc(distance,frm.Time,lp(frm.Data));
      hold on; plot(distance,sf,'--'); plot(distance,bttime_n,'--');
  end;
 

%% Correction for Elevation
sf_new = zeros(size(sf));
bt_new = zeros(size(bttime_n));
frm.Elevation_new = zeros(size(frm.Elevation));
frm.sf_elev_new = zeros(size(frm.Elevation));
frm.bt_elev_new = zeros(size(frm.Elevation));

% 1)Remove data before zero time
negative_bins = frm.Time < 0;
frm.Time_new = frm.Time(~negative_bins);
frm.Data_new = frm.Data(~negative_bins,:);

% 2)Create elevation axis to interpolate to
[max_elev,max_elev_idx] = max(frm.Elevation);
min_elev = min(frm.Elevation - sf*c/2 - (frm.Time_new(end)-sf)*c/2/sqrt(er_ice));
dt = frm.Time(2)-frm.Time(1);
dr = dt * c/2 / sqrt(er_ice);
dt_air = dr/(c/2);
elev_axis = max_elev:-dr:min_elev;
new_time = zeros(length(elev_axis),length(frm.Elevation));

% 3)Zero pad data to create space for interpolated data
zero_pad_len = length(elev_axis) - length(frm.Time_new);
frm.Data_new = cat(1,frm.Data_new,zeros(zero_pad_len,size(frm.Data_new,2)));

% 4)Determine the corrections to apply to elevation and layers
dRange = max_elev - frm.Elevation;
dBins = round(dRange / (c/2) / dt);
dtime = dRange/(c/2);

for rline = 1:size(frm.Data_new,2)
    % Determine elevation bins before surface
    sf_elev = frm.Elevation(rline) - sf(rline) * c/2;
    time0 = -(max_elev - frm.Elevation(rline))/(c/2);
    last_air_idx = find(elev_axis > sf_elev,1,'last');
    new_time_tmp = (time0 + dt_air*(0:last_air_idx-1)).';
    if last_air_idx < length(elev_axis)
        % Determine elevation bins after surface
        dt_ice = dr/(c/2/sqrt(er_ice));
        first_ice_idx = last_air_idx + 1;
        time0 = sf(rline) + (sf_elev - elev_axis(first_ice_idx))/(c/2/sqrt(er_ice));
        new_time(:,rline) = cat(1,new_time_tmp, (time0 + dt_ice*(0:length(elev_axis)-length(new_time_tmp)-1)).');
    end
    frm.Data_new(:,rline) = interp1(frm.Time_new, frm.Data_new(1:length(frm.Time_new),rline), new_time(:,rline), 'linear',0);
    frm.Elevation_new(rline) = frm.Elevation(rline) + dRange(rline);
    sf_new(rline) = sf(rline) + dtime(rline);
    bt_new(rline) = bttime_n(rline) + dtime(rline);
    frm.sf_elev_new(rline) = frm.Elevation_new(rline) - sf_new(rline)*c/2;
    frm.bt_elev_new(rline) = frm.sf_elev_new(rline) - (bt_new(rline)-sf_new(rline))*c/2/sqrt(er_ice);
end

if debug_flag==1
 fh = figure(2);
 figure(fh);imagesc([],elev_axis,10*log10(abs(frm.Data_new).^2));title('Elevation Correction Data')
 ax = gca;
 ax.YDir = 'normal';
 hold on;plot(frm.sf_elev_new,'--');plot(frm.bt_elev_new,'--');
 bt_slope = diff(frm.bt_elev_new)./diff(distance);
end

%% Ice Slope Correction
p=4.99; %Prewindowed gate 
MeanDepth=(-nanmean(sf)+nanmean(bttime_n))*c/(2*sqrt(er_ice));    %Mean Ice Depth 
Nx_int_dist=2*sqrt(MeanDepth*p/sqrt(er_ice));                    %Incoherent Integration distance
Nx_int = floor(Nx_int_dist/(distance(10)-distance(9)));                %No of lines integrated
slopeval=zeros(1,size(frm.Data_new,1));
Nx0 = size(frm.Data_new,2);
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
    
    if isinf(nanmean(frm.bt_elev_new(idx1:idx2))) || isnan(nanmean(frm.bt_elev_new(idx1:idx2))) % If No Ice Bottom, skip
         continue;
    end
    
    p = polyfit(distance(idx1:idx2),frm.bt_elev_new(idx1:idx2),1);    % Polygonal fitting
    slopeval(idx1:idx2) = polyval(p,distance(idx1:idx2));           % Polygonal values after fitting
    base = distance(idx2)-distance(idx1);
    perpendicular = slopeval(idx2)-slopeval(idx1);
    hypotenuse = sqrt(base^2+perpendicular^2);
    angle1(rline)=asin(perpendicular/hypotenuse)*180/pi;
    slopeerror = slopeval(idx1:idx2)-slopeval(idx1);         % Error betn original and fitting line
    dtime = 2*slopeerror/c/sqrt(er_ice);
    
    for idx = idx1:idx2
        frm.Data_new(:,idx) = interp1(elev_axis, frm.Data_new(:,idx), elev_axis + slopeerror(idx-idx1 +1), 'linear',0);
        frm.Elevation_new(idx) = frm.Elevation_new(idx) - slopeerror(idx-idx1 +1);
        frm.sf_elev_new(idx) = frm.sf_elev_new(idx) - slopeerror(idx-idx1 +1);
        frm.bt_elev_new(idx) = frm.bt_elev_new(idx) - slopeerror(idx-idx1 +1);  
        sf_new(idx) = sf_new(idx) + dtime(idx-idx1 +1);
        bt_new(idx) = bt_new(idx) + dtime(idx-idx1 +1);
    end
   
end

if debug_flag==1
 fh = figure(3);
 figure(fh);imagesc([],elev_axis,10*log10(abs(frm.Data_new).^2));title('Slope Correction Data')
 ax = gca;
 ax.YDir = 'normal';
 hold on;plot(frm.sf_elev_new,'--');plot(frm.bt_elev_new,'--');
end

%%
 % Truncate data around ice bottom within bt.range_bins
  bt_range_bins =param.analysis.bt.range_bins;
  bt.peakval = NaN*ones(1,size(frm.Data_new,2));
  bt.peakidx = NaN*ones(1,size(frm.Data_new,2));
  bt.waveform = NaN*ones(length(bt_range_bins),size(frm.Data_new,2));
  bt.inc_wf_ave=NaN*ones(size(bt.waveform,1),Nx);
  for rline = 1:size(frm.Data_new,2)
    if ~isnan(frm.bt_elev_new(rline)) & ~isinf(frm.bt_elev_new(rline)) 
     % bt.idx(rline) = find(elev_axis<=frm.bt_elev_new(rline),1,'first');
     bt.peakidx(rline) = round(interp1(elev_axis,[1:length(elev_axis)],frm.bt_elev_new(rline)));
     bt.peakval(rline) = frm.Data_new(bt.peakidx(rline),rline);
      first_idx = bt.peakidx(rline) + bt_range_bins(1); 
      last_idx = bt.peakidx(rline) + bt_range_bins(end);
      if first_idx < 1 | last_idx>size(frm.Data_new,1)
        bt.peakidx(rline) = NaN;
        bt.peakval(rline) = NaN;
        continue
      end
      lower=bt.peakidx(rline)+bt_range_bins(1);
      upper=bt.peakidx(rline)+bt_range_bins(end);
      if upper >size(frm.Data_new,1)
          upper=size(frm.Data_new,1);
      end
      bt.waveform(1:51+upper-bt.peakidx(rline),rline) = frm.Data_new(lower:upper, rline);
    else
      continue
    end
  end
 
 if debug_flag==1
    figure(4);imagesc(lp(frm.Data_new));
    hold on;plot(bt.peakidx,'--');
    figure(5);plot(lp(bt.waveform));
    figure(6);imagesc(lp(bt.waveform));
 end
 

%% Coherene Index Calculation with Slope error Correction


bt.coherenceindex.Latitude_mean=zeros(1,Nx);
bt.coherenceindex.Longitude_mean=zeros(1,Nx);
bt.coherenceindex.GPS_time_ave=zeros(1,Nx);
square_int = zeros(size(frm.Data_new,1),Nx); % incoherent integration, take square of abs first, then sum 
square_int_dB=zeros(size(frm.Data_new,1),Nx);
int_square = zeros(size(frm.Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
coh_index_slopecorr = NaN*ones(1,Nx);
Abruptiveindex=NaN*ones(1,Nx);
risingedge=NaN*ones(1,Nx);
fallingedge=NaN*ones(1,Nx);
bt.avepeakindex=NaN*ones(1,Nx);
Padj=NaN*ones(1,Nx);

for rline = 1:Nx
    idx1 = (rline-1)*Nx_int + 1;
    idx2 = rline*Nx_int;
    if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
        idx2 = Nx0;
    else
        idx2 = min(idx2,Nx0);
    end
    
    bt.coherenceindex.Latitude_mean(rline)=mean(frm.Latitude(idx1:idx2));  %Mean Latitude
    bt.coherenceindex.Longitude_mean(rline)=mean(frm.Longitude(idx1:idx2));
    bt.coherenceindex.GPS_time_ave(rline)=mean(frm.GPS_time(idx1:idx2));
    
    square_int(:,rline) = mean(abs(frm.Data_new(:,idx1:idx2)).^2,2);
    int_square(:,rline) = abs(mean(frm.Data_new(:,idx1:idx2),2)).^2;
    square_int_dB(:,rline) = 10*log10(square_int(:,rline));
    
    meanbt=nanmean(frm.bt_elev_new(idx1:idx2)); 
    if meanbt==0 || isnan(meanbt)
        continue;                %skip if no ice bottom
    end
   
    meansf_og=nanmean(sf(idx1:idx2));   %If bottm=surface then skip
    meanbt_og=nanmean(bttime(idx1:idx2));
    if meanbt_og==meansf_og
        continue;
    end
    %bt_idx_m = find(elev_axis<=meanbt,1,'first');
    bt_idx_m = round(interp1(elev_axis,[1:length(elev_axis)],meanbt));

    if isempty(bt_idx_m)
      continue;      %skip if no ice bottom 
    end
     
    b1=bt_idx_m+param.analysis.bt.peak_bins(1);               % Peak search index peak bins
    b2=bt_idx_m+param.analysis.bt.peak_bins(2);
    if b2>size( frm.Data_new,1)
        b2=size( frm.Data_new,1);
    end
     
    [bt_val,bt_idx] = max(square_int(b1:b2,rline));   %Peak Index and Value
    bt_idx = bt_idx_m+param.analysis.bt.peak_bins(1)+bt_idx-1;
    bt_pwr = 10*log10(bt_val);
    noise_bin1=bt_idx+70;
    noise_bin2=bt_idx+70+param.analysis.bt.noise_bins;
    if noise_bin1>size( frm.Data_new,1) || noise_bin2>size( frm.Data_new,1)
        noise_bin1=size( frm.Data_new,1)-param.analysis.bt.noise_bins;
        noise_bin2=size( frm.Data_new,1);
    end
   
   noise = mean(square_int_dB(noise_bin1:noise_bin2,rline));
   if bt_pwr-noise<3
      % noise_bin1=size( frm.Data_new,1)-100;
       %noise_bin2=size( frm.Data_new,1)-70;
       continue;
   end
   %noise = mean(square_int_dB(noise_bin1:noise_bin2,rline));
   risingedge(rline) = bt_idx-1;
   SNR = bt_pwr-noise;
    
    % find risingedge(rline): risingedge and fallingedge(rline):fallingedge for integration in range
   while square_int_dB(risingedge(rline),rline)-noise > 0.05*SNR  && risingedge(rline)>2
        risingedge(rline) = risingedge(rline) - 1;
   end
   if bt_idx-risingedge(rline)>param.analysis.bt.max_index_range
        risingedge(rline)=bt_idx-param.analysis.bt.max_index_range;  %Max 50 bins away for risingedge
   end
   
   if bt_idx-risingedge(rline)<param.analysis.bt.guard_bins
       risingedge(rline)=bt_idx-param.analysis.bt.guard_bins;  %Guard band of 3
   end
   
   [r,tmp_idx]=(min(square_int_dB(risingedge(rline):bt_idx-param.analysis.bt.guard_bins,rline)));
   risingedge=risingedge+tmp_idx-1;
   
   fallingedge(rline) = bt_idx+1;  
   if fallingedge(rline)>size( frm.Data_new,1)
       fallingedge(rline)=size( frm.Data_new,1);
   end
   while square_int_dB(fallingedge(rline),rline)-noise > 0.05*SNR && fallingedge(rline)<size(frm.Data_new,1) 
      fallingedge(rline) = fallingedge(rline) + 1;                                 
   end
   if fallingedge(rline)-bt_idx>param.analysis.bt.max_index_range
        fallingedge(rline)=bt_idx+param.analysis.bt.max_index_range;  %Max 50 bins away for falling edge
   end
   if fallingedge(rline)-bt_idx<param.analysis.bt.guard_bins
        fallingedge(rline)=bt_idx+param.analysis.bt.guard_bins;
   end
   if fallingedge(rline)>size( frm.Data_new,1)
       fallingedge(rline)=size( frm.Data_new,1);
   end
   coh_index_slopecorr(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline)); %Coherence Index
     
   Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));  %Power summed from risingedge(rline) to fallingedge(rline)
   Abruptiveindex(rline)=bt_val/Imeanx;    %Abruptive Index
    
   bt.avepeakindex(rline)=bt_idx;           %peakindex
   bt.inc_wf_ave(:,rline) = square_int_dB(bt_idx+bt_range_bins);  %Average Waveform
   
    B=2.3*3000/(depth+2000);
    Padj(rline)=(lp(Imeanx*depth.^2)+B*depth/100-2); 
end

if debug_flag==1
    hold on; figure(8); hold on; plot(coh_index_slopecorr,'g','Displayname','SlopeCorrected'); title('Coherence Index')
    legend('show');
    figure(7);plot(coh_index_slopecorr),title('Coherence Index after Slope Correction')
end

%%


bt.coherenceindex.GPS_time=frm.GPS_time;
bt.coherenceindex.Latitude=frm.Latitude;
bt.coherenceindex.Longitude=frm.Longitude;
bt.coherenceindex.Elevation_new=frm.Elevation_new;
bt.coherenceindex.Surface_new=frm.sf_elev_new;
bt.coherenceindex.SurfaceTime_new=sf_new;
bt.coherenceindex.Bottom_new=frm.bt_elev_new;
bt.coherenceindex.BottomTime_new=bt_new;
bt.coherenceindex.Time_new=frm.Time_new';
bt.coherenceindex.frm_id={};
bt.coherenceindex.value=coh_index_slopecorr;
bt.coherenceindex.Abruptiveindex=Abruptiveindex;
bt.coherenecindex.AdjustedIntensity=Padj;

end

% Save
out_fn = fullfile(ct_filename_out(param, ...
  param.analysis.out_path,'CSARP_basal_condition'), ...
  sprintf('IceBedCoherenceIndex_%s_%03d.mat',param.day_seg,param.load.frm));

[out_fn_dir] = fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end

param_analysis = param;
fprintf('  Saving outputs %s\n', out_fn);
save(save_format, out_fn, 'bt','param_analysis');

success = true;

return;
