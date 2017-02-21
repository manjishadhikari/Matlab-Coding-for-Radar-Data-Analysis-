
close all;
load_flag=0;
if load_flag
     frm= load('/cresis/snfs1/dataproducts/ct_data/rds/2011_Antarctica_TO/CSARP_ant/20111201_04/Data_20111201_04_002.mat');
    layer = load('/cresis/snfs1/dataproducts/public/data/rds/2011_Antarctica_TO/CSARP_layerData/20111201_04/Data_20111201_04_002.mat');
end
distance=geodetic_to_along_track(frm.Latitude,frm.Longitude,frm.Elevation);
  
  fc = 195e6;
fs = 1.1111e8;  %Sampling Frequency of the radar
c = 3e8;        %Speed of light
er_ice = 3.15;  %Permittivity of ice

 frm.Bottom = layer.layerData{2}.value{2}.data;

  
    if length(frm.Bottom) ~= size(frm.Data,2) 
      if all(isnan(frm.Bottom))
        frm.Bottom = NaN*ones(size(frm.Data,1));
      else
        nan_idxs = find(isnan(frm.Bottom));
        layer.GPS_time(nan_idxs) = [];
        frm.Bottom(nan_idxs) = [];
        frm.Bottom = interp1(layer.GPS_time,frm.Bottom,frm.GPS_time,'linear','extrap');
      end
    end

  
  frm.Bottom(frm.Bottom>frm.Time(end) | frm.Bottom<frm.Time(1)) = NaN; % remove pick errors exist
  debug_flag=1;
  if debug_flag
      figure(1); imagesc(frm.GPS_time,frm.Time,lp(frm.Data));
      hold on; plot (frm.GPS_time,frm.Bottom,'--');
  end;

  % Truncate data around ice bottom within bt.range_bins
  bt_range_bins =[-50:100];
  bt.val = NaN*ones(1,size(frm.Data,2));
  bt.idx = NaN*ones(1,size(frm.Data,2));
  bt.waveform = NaN*ones(length(bt_range_bins),size(frm.Data,2));
  for rline = 1:size(frm.Data,2)
    if ~isnan(frm.Bottom(rline)) & ~isinf(frm.Bottom(rline)) 
      bt.idx(rline) = round(interp1(frm.Time,[1:length(frm.Time)],frm.Bottom(rline)));
      bt.val(rline) = frm.Data(bt.idx(rline),rline);
      first_idx = bt.idx(rline) + bt_range_bins(1); 
      last_idx = bt.idx(rline) + bt_range_bins(end);
      if first_idx < 1 | last_idx>size(frm.Data,1)
        first_idx=NaN;
        last_idx=NaN;
        continue;
      end
      bt.waveform(:,rline) = frm.Data(bt.idx(rline) + bt_range_bins, rline);
    else
      continue
    end
  end

  if 1
    figure(3);imagesc(lp(frm.Data));
    hold on;plot(bt.idx,'--');
    figure(4);plot(lp(bt.waveform));
  end

  % average N ice bed echo waveforms, for specular area, we assume the
  % averaged waveform would have high SNR and narrower 6db pulse widith
  % becausee of coherence.
  
  rising_edge_bins = -50 : -3;
  falling_edge_bins = 3 : 100;
  specular_threshold = 40;
  specular_pulse_width = 4;
  specular_rising_edge_SL =3;
  specular_falling_edge_SL = 20;
  param.analysis.bt.specular_SL_guard_bins=3;
  param.analysis.bt.peak_bins=[-5:5];
  N =750;
  blk_N = floor(size(frm.Data,2)/N);
  bt.GPS_time_ave = NaN*ones(1,blk_N);
  bt.SNR_ave = NaN*ones(1,blk_N);
  bt.pulse_width_ave = NaN*ones(1,blk_N);
  bt.rising_edge_SL = NaN*ones(1,blk_N);
  bt.falling_edge_SL = NaN*ones(1,blk_N);
  bt.wf_ave = NaN*ones(length(bt_range_bins),blk_N);
  bt.specular.GPS_time =  [];
  bt.specular.Longitude = [];
  bt.specular.Latitude = [];
  bt.specular.Elevation = [];
  bt.specular.SNR = [];
  bt.specular.pulse_width = [];
  bt.specular.rising_edge_SL = [];
  bt.specular.falling_edge_SL = [];
  bt.specular.Bottom = [];
  bt.specular.Surface = [];
  bt.specular.depth = [];
  bt.specular.waveform = [];
 
  square_int = zeros(size( bt.wf_ave,1), blk_N); % incoherent integration, take square of abs first, then sum 
int_square = zeros(size(bt.wf_ave,1), blk_N); % coherent integration, sum complex data first, then take square of abs
coh_index_slopecorr =NaN*ones(1, blk_N);
 bt.fallingedge=zeros(size(1,blk_N));
 bt.fallingedge=zeros(size(1,blk_N));
 bt.meanidx=zeros(1,blk_N);
  for blk_idx =1:blk_N
     idx1=(blk_idx-1)*N+1 ;
    idx2=blk_idx*N;
      bt.wf_ave(:,blk_idx) = mean(bt.waveform(:,(blk_idx-1)*N +1:blk_idx*N),2);
    [bt_peak,bt_peak_idx] = max(bt.wf_ave(1-bt_range_bins(1)+ param.analysis.bt.peak_bins,blk_idx));
    bt_peak_idx = bt_peak_idx - bt_range_bins(1) -5;
    bt.SNR_ave(blk_idx) = lp(bt_peak/mean(bt.wf_ave(end-30:end,blk_idx)));
    bt.GPS_time_ave(blk_idx) = frm.GPS_time((blk_idx-1)*N +round(0.5*N));
    rising_idx = bt_peak_idx-1;
    while rising_idx >= 2 && lp(bt.wf_ave(rising_idx,blk_idx)) > lp(bt_peak)-6
      rising_idx = rising_idx - 1;
    end
    falling_idx = bt_peak_idx+1;
    while falling_idx < length(bt.wf_ave(:,blk_idx)) && lp(bt.wf_ave(falling_idx)) > lp(bt_peak)-6
      falling_idx = falling_idx + 1;
    end
    bt.pulse_width_ave(blk_idx) = falling_idx -rising_idx;
    bt.rising_edge_SL(blk_idx) = lp(max(bt.wf_ave(1:max(1,bt_peak_idx-param.analysis.bt.specular_SL_guard_bins),blk_idx))) - lp(bt_peak);    
    bt.falling_edge_SL(blk_idx) = lp(max(bt.wf_ave(min(bt_peak_idx+param.analysis.bt.specular_SL_guard_bins,size(bt.wf_ave,1)):end,blk_idx))) - lp(bt_peak);
    bt.fallingedge(blk_idx)=falling_idx;
    bt.meanidx(blk_idx)=bt_peak_idx;
    
  end
 
  
  if 1
    figure(5);plot(lp(bt.wf_ave));
    figure(6);subplot(2,1,1);plot(bt.SNR_ave);subplot(2,1,2);plot(bt.pulse_width_ave);
    figure(7);imagesc(frm.GPS_time,frm.Time,lp(frm.Data));
    specular_bt = interp1(frm.GPS_time,frm.Bottom,bt.specular.GPS_time);
    hold on;plot(bt.specular.GPS_time,specular_bt,'o');
    figure(8);subplot(5,1,1);plot(bt.specular.SNR);
    subplot(5,1,2);plot(bt.specular.pulse_width);
    subplot(5,1,3);plot(bt.specular.rising_edge_SL);
    subplot(5,1,4);plot(bt.specular.falling_edge_SL);
    subplot(5,1,5);plot(bt.specular.depth);
    figure(9);plot(lp(bt.specular.waveform))
  end

  for blk_idx = 1:Nx
      
      idx1=(blk_idx-1)*Nx_int +1;
      idx2= blk_idx*Nx_int;   
      if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
        idx2 = Nx0;
      else
        idx2 = min(idx2,Nx0);
      end
      meanbt=nansum(bt.waveform(isfinite(bt.waveform(idx1:idx2)))); 
      count=nansum(isfinite(bt.waveform(idx1:idx2)));
      bt.wf_ave=meanbt/count;
   
  end
  
    bt_waveforms=bt.waveform(idx1:idx2);
     sum_bt_waveform=nansum(bt_waveforms(isfinite(bt_waveforms))); 
     count=nansum(isfinite(bt_waveforms));
     bt.wf_ave(:,rline)=sum_bt_waveform/count;
     
      bt.wf_ave(:,rline) = nanmean(A.Data_new(bt_idx-50:bt_idx+100,idx1:idx2),2)