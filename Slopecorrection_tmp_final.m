
%% Determination of coherence Index 

tic
fc = 195e6;
fs = 1.1111e8;  %Sampling Frequency of the radar
c = 3e8;        %Speed of light
er_ice = 3.15;  %Permittivity of ice
p=4.99;   %Pre windowed radar pulse halfwidth in air     

load_flag = 0;

if load_flag
     A = load('/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_CSARP_manjish/20120404_01/Data_20120404_01_046.mat');
    L = load('/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_post/CSARP_layerData/20120404_01/Data_20120404_01_046.mat');
end
toc

distance1=geodetic_to_along_track(A.Latitude,A.Longitude,A.Elevation); %Horizontal Distance
sf1 = interp1(L.GPS_time,L.layerData{1}.value{2}.data,A.GPS_time,'linear','extrap'); %Two way prop time to surface
bttime1= interp1(L.GPS_time, L.layerData{2}.value{2}.data,A.GPS_time,'linear','extrap');
bttime1(bttime1>A.Time(end) | bttime1<A.Time(1)) = NaN; % remove pick errors exist


figure(40),imagesc(distance1,A.Time,10*log10(abs(A.Data)));
hold on;plot(distance1,sf1,'--');plot(distance1,bttime1,'--');

Nx_int=floor(4/distance1(2));
 Nx0 = size(A.Data,2);
   Nx = floor(Nx0/Nx_int); 
 Nx_mod = mod(Nx0,Nx_int);
 if Nx_mod>= Nx_int/2;
      Nx = Nx + 1;
 end
 A.Dataweiner=NaN*ones(size(A.Data,1),Nx);
 distance=NaN*ones(1,Nx);
 bttime=NaN*ones(1,Nx);
 sf=NaN*ones(1,Nx);
  A.Elevation1=NaN*ones(1,Nx);
for rline =1:Nx
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        A.Dataweiner(:,rline)=mean(A.Data(:,idx1:idx2),2);
        distance(rline)=mean(distance1(:,idx1:idx2),2);
        bttime(rline)=mean(bttime1(:,idx1:idx2),2);
        sf(rline)=mean(sf1(:,idx1:idx2),2);
         A.Elevation1(:,rline)=mean(A.Elevation(:,idx1:idx2),2);
end

%bttime=[bttime(1:6000),bttime(1:end-6000)]; 
%bttime=bttime-0.015E-05;
figure(43),imagesc(distance,A.Time,10*log10(abs(A.Data)));
hold on;plot(distance,sf,'--');plot(distance,bttime,'--');
%A.Dataweiner=wiener2(A.Dataweiner,[3,10]);
% f= fspecial('gaussian',3,100);
% A.Dataweiner= imfilter(A.Data,f,'replicate');

%A.Dataweiner=mean(A.Data(:,1:8:end),1);
%distance=distance(:,1:8:end);
%A.Dataweiner=A.Data;

 debug_flag=0;
  if debug_flag
      figure(10); imagesc(A.GPS_time,A.Time,lp(A.Dataweiner));
      hold on; plot (A.GPS_time,bttime,'--');
  end; 
  power=lp(abs(A.Dataweiner));
  
%bttime_n=bttime;  
bttime=[bttime(1:80),bttime(1:end-80)];  
% 
bttime_n=NaN*ones(1,size(A.Dataweiner,2));
for rline=1:size(A.Dataweiner,2)
   if isnan(bttime1(rline)) || isinf(bttime1(rline))
       bttime(rline)=NaN;
       continue;
   end
    bt_idx_m = find(A.Time>bttime(rline),1); 
    if bt_idx_m-15>0 & bt_idx_m+15<size(A.Dataweiner,1)
    [pwr,bt_idx]=max(power(bt_idx_m-15:bt_idx_m+15,rline));
    dt=A.Time(2)-A.Time(1);
    bttime_n(rline)=(bttime(rline)-15*dt)+(bt_idx-1)*dt;
    else continue;
    end
    if bttime_n(rline)<sf(rline)
        bttime_n(rline)=sf(rline);
    end
end

%bttime_n=bttime_n+0.001*10^-5;
%bttime=sgolayfilt(bttime,1,251);
% sf=sf(:,1:8:end);
% bttime_n=bttime_n(:,1:8:end);

figure(1),imagesc(distance,A.Time,10*log10(abs(A.Dataweiner)));
hold on;plot(distance,sf,'--');plot(distance,bttime_n,'--');

%figure(1),imagesc([1:size(distance,2)],A.Time,10*log10(abs(A.Dataweiner)));
%hold on;plot([1:size(distance,2)],sf,'--');plot([1:size(distance,2)],bttime,'--');


%% Coherene Index Calculation with original Data
sf(sf==Inf)=NaN;
bttime_n(bttime_n==Inf)=NaN;
MeanDepth_og=(nanmean(bttime)-nanmean(sf))*c/(2*sqrt(er_ice));  %Original Mean Depth
MeanDepth=(nanmean(bttime_n)-nanmean(sf))*c/(2*sqrt(er_ice)); 
 %Mean Ice Depth 

if MeanDepth_og~=0 & MeanDepth_og>150;
    Nx_int_dist=2*sqrt(MeanDepth*p/sqrt(er_ice));                    %Incoherent Integration distance
    Nx_int = floor(Nx_int_dist/(distance(10)-distance(9)));                %No of lines integrated
    Nx0 = size(A.Dataweiner,2);
    Nx = floor(Nx0/Nx_int);                                         %Total lines after integration
    Nx_mod = mod(Nx0,Nx_int);
    if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
    end
    square_int = zeros(size(A.Dataweiner,1),Nx); % incoherent integration, take square of abs first, then sum (phase info lost)
    int_square = zeros(size(A.Dataweiner,1),Nx); % coherent integration, sum complex data first, then take square of abs (phase info remains)
    coh_index_ogdata = zeros(1,Nx);
    Abruptiveindex=zeros(1,Nx);
    risingedge=zeros(1,Nx);
    fallingedge=zeros(1,Nx);
    Padj=zeros(1,Nx);
    for rline =1:Nx
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        square_int(:,rline) = mean(abs(A.Dataweiner(:,idx1:idx2)).^2,2);
        square_int_dB = 10*log10(square_int(:,rline));
        int_square(:,rline) = (abs(mean(A.Dataweiner(:,idx1:idx2),2))).^2;
        
        meanbt=nansum(bttime_n(isfinite(bttime_n(idx1:idx2))));
        count=nansum(isfinite(bttime_n(idx1:idx2)));
        if meanbt==0 || count==0
            continue;
        end
        meanbt=meanbt/count;
        bt_idx_m = find(A.Time>meanbt,1,'first');
        
        if isempty(bt_idx_m)
            bt_idx_m=round(size(A.Dataweiner,1)/2);
        end
        b1=bt_idx_m-5;
        b2=bt_idx_m+5;
        if b2>size(A.Dataweiner,1)
            b2=size(A.Dataweiner,1);
        end
        
        [bt_val,bt_idx] = max(square_int(b1:b2,rline));
        bt_idx = bt_idx_m-5+bt_idx-1;                 %Peak Index
        bt_pwr = 10*log10(bt_val);                     %Peak Value
        square_int_dB = 10*log10(square_int(:,rline));
        noise_bin1=bt_idx+70;
        noise_bin2=bt_idx+100;
        if noise_bin1>size(A.Dataweiner,1)
            noise_bin1=size(A.Dataweiner,1)-30;
            noise_bin2=size(A.Dataweiner,1);
        end
        
        if noise_bin2>size(A.Dataweiner,1)
            noise_bin1=size(A.Dataweiner,1)-30;
            noise_bin2=size(A.Dataweiner,1);
        end
        noise =mean(square_int_dB(noise_bin1:noise_bin2));
        risingedge(rline) = bt_idx-1;
        SNR = bt_pwr-noise;
        
        % find risingedge(rline) and fallingedge(rline) for integration in range
        while square_int_dB(risingedge(rline))-noise > 0.05*SNR && risingedge(rline)>2
            risingedge(rline) = risingedge(rline) - 1;
        end
        if risingedge(rline)==0
            risingedge(rline)=find(square_int_dB==min(square_int_dB(b1:bt_idx)),1,'first');
        end
        
        fallingedge(rline) = bt_idx+1;
        if fallingedge(rline)>size(A.Dataweiner,1)
            fallingedge(rline)=size(A.Dataweiner,1) ;
        end
        while square_int_dB(fallingedge(rline))-noise > 0.05*SNR &&  fallingedge(rline)<size(A.Dataweiner,1)      %Depth Bins risingedge(rline) and fallingedge(rline)
            fallingedge(rline) = fallingedge(rline) + 1;
        end
        if fallingedge(rline)-bt_idx>50        %Max Allowed bin
            fallingedge(rline)=bt_idx+50;
        end
        if fallingedge(rline)-bt_idx<3           %Guard bin
            fallingedge(rline)=bt_idx+3;
        end
        D1=risingedge(rline);
        D2=fallingedge(rline);
        
        Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));
        Abruptiveindex(rline)=bt_val/Imeanx;
        coh_index_ogdata(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline));
   
    
    end
    distancenx=max(distance).*Nx_int*0.5/(size(distance,2)):max(distance).*Nx_int/size(distance,2):max(distance).*Nx_int/size(distance,2)*Nx;
    figure(4);plot(distancenx,coh_index_ogdata,'Displayname','Original Data');
    title('Original Data Coherence Index');   %Coherence Index without Correction
    figure(40);plot(distancenx,Abruptiveindex);
   
    %% Correction for Elevation
    sf_new = zeros(size(sf));
    bt_new = zeros(size(bttime_n));
    A.Elevation1_new = zeros(size(A.Elevation1));
    A.sf_elev_new = zeros(size(A.Elevation1));
    A.bt_elev_new = zeros(size(A.Elevation1));
    
    % 1)Remove data before zero time
    negative_bins = A.Time < 0;
    A.Time_new = A.Time(~negative_bins);
    A.Data_new = A.Dataweiner(~negative_bins,:);
    
    % 2)Create elevation axis to interpolate to
    [max_elev,max_elev_idx] = max(A.Elevation1);
    min_elev = min(A.Elevation1 - sf*c/2 - (A.Time_new(end)-sf)*c/2/sqrt(er_ice));
    dt = A.Time(2)-A.Time(1);
    dr = dt * c/2 / sqrt(er_ice);
    dt_air = dr/(c/2);
    elev_axis = max_elev:-dr:min_elev;
    new_time = zeros(length(elev_axis),length(A.Elevation1));
    
    % 3)Zero pad data to create space for interpolated data
    zero_pad_len = length(elev_axis) - length(A.Time_new);
    A.Data_new = cat(1,A.Data_new,zeros(zero_pad_len,size(A.Data_new,2)));
    
    % 4)Determine the corrections to apply to elevation and layers
    dRange = max_elev - A.Elevation1;
    dBins = round(dRange / (c/2) / dt);
    dtime = dRange/(c/2);
    
    for rline = 1:size(A.Data_new,2)
        % Determine elevation bins before surface
        sf_elev = A.Elevation1(rline) - sf(rline) * c/2;
        time0 = -(max_elev - A.Elevation1(rline))/(c/2);
        last_air_idx = find(elev_axis > sf_elev,1,'last');
        new_time_tmp = (time0 + dt_air*(0:last_air_idx-1)).';
        if last_air_idx < length(elev_axis)
            % Determine elevation bins after surface
            dt_ice = dr/(c/2/sqrt(er_ice));
            first_ice_idx = last_air_idx + 1;
            time0 = sf(rline) + (sf_elev - elev_axis(first_ice_idx))/(c/2/sqrt(er_ice));
            new_time(:,rline) = cat(1,new_time_tmp, (time0 + dt_ice*(0:length(elev_axis)-length(new_time_tmp)-1)).');
        end
        A.Data_new(:,rline) = interp1(A.Time_new, A.Data_new(1:length(A.Time_new),rline), new_time(:,rline), 'linear',0);
        A.Elevation1_new(rline) = A.Elevation1(rline) + dRange(rline);
        sf_new(rline) = sf(rline) + dtime(rline);
        bt_new(rline) = bttime_n(rline) + dtime(rline);
        A.sf_elev_new(rline) = A.Elevation1_new(rline) - sf_new(rline)*c/2;
        A.bt_elev_new(rline) = A.sf_elev_new(rline) - (bt_new(rline)-sf_new(rline))*c/2/sqrt(er_ice);
          
       
        
    end
    
    fh = figure(2);
    figure(fh);imagesc(distance,elev_axis,10*log10(abs(A.Data_new).^2));title('Elevation Correction Data')
    ax = gca;
    ax.YDir = 'normal';
    hold on;plot(distance,A.sf_elev_new,'--');plot(distance,A.bt_elev_new,'--');
    bt_slope = diff(A.bt_elev_new)./diff(distance);
   
    
    
    
    %% Coherene Index Calculation with elevation Correction
    
    square_int = zeros(size( A.Data_new,1),Nx); % incoherent integration, take square of abs first, then sum
    int_square = zeros(size(A.Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
    coh_index_elevcorr = zeros(1,Nx);
    for rline =1:Nx
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        square_int(:,rline) = mean(abs(A.Data_new(:,idx1:idx2)).^2,2);
        int_square(:,rline) = abs(mean(A.Data_new(:,idx1:idx2),2)).^2;
        square_int_dB=lp(square_int(:,rline));
        
        meanbt=nansum(A.bt_elev_new(idx1:idx2));
        count=nansum(isfinite(A.bt_elev_new(idx1:idx2)));
        meanbt=meanbt/count;
        if meanbt==0 || count==0 ||isnan(meanbt)
            continue;
        end
        
        bt_idx_m = find(elev_axis<=meanbt,1,'first');
        
        if isempty(bt_idx_m)
            bt_idx_m=round(size(A.Data_new,1)/2);      %looking for ice bottom
        end
        b1=bt_idx_m-5;
        b2=bt_idx_m+5;
        if b2>size(A.Data_new,1)
            b2=size(A.Data_new,1);
        end
        
        [bt_val,bt_idx]=max(square_int(b1:b2,rline));
        bt_idx = bt_idx_m-5+bt_idx-1;
        bt_pwr = 10*log10(bt_val);
        noise_bin1=bt_idx+70;
        noise_bin2=bt_idx+100;
        if noise_bin1>size(A.Data_new,1)
            noise_bin1=size(A.Data_new,1)-30;
            noise_bin2=size(A.Data_new,1);
        end
        if noise_bin2>size(A.Data_new,1)
            noise_bin2=size(A.Data_new,1);
        end
        noise = mean(square_int_dB(noise_bin1:noise_bin2));
        
        risingedge(rline) = bt_idx-1;
        SNR = bt_pwr-noise;
        
        % find risingedge(rline) and fallingedge(rline) for integration in range
        while square_int_dB(risingedge(rline))-noise > 0.05*SNR & risingedge(rline)>2
            risingedge(rline) = risingedge(rline) - 1;
        end
        if risingedge(rline)==0|square_int_dB(risingedge(rline))==-Inf|square_int_dB(risingedge(rline))==Inf
            risingedge(rline)=find(square_int_dB==min(square_int_dB(b1:bt_idx)),1,'first');
        end
        fallingedge(rline) = bt_idx+1;
        if fallingedge(rline)>size(A.Data_new,1)
            fallingedge(rline)=size(A.Data_new,1) ;
        end
        while square_int_dB(fallingedge(rline))-noise > 0.05*SNR && fallingedge(rline)<size(A.Data_new,1)      %Depth Bins risingedge(rline) and fallingedge(rline)
            fallingedge(rline) = fallingedge(rline) + 1;
        end
        Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));
        Abruptiveindex(rline)=bt_val/Imeanx;
        coh_index_elevcorr(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline));
    end
    hold on; figure(4); hold on; plot(distancenx,coh_index_elevcorr,'r','Displayname','Elevation Corrected');
    figure(41); hold on; plot(distancenx,Abruptiveindex,'r','Displayname','Elevation Corrected'); grid;
    
    %% Ice Slope Correction
    
    Nx0 = size(A.Data_new,2);
    Nx = floor(Nx0/Nx_int);
    Nx_mod = mod(Nx0,Nx_int);
    if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
    end
    slopeerror=zeros(1,size(A.Data_new,1));
    slopeval=zeros(1,size(A.Data_new,1));
    angle=zeros(1,Nx);
    
    for rline =1:Nx
        
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        
        if isinf(mean(A.bt_elev_new(idx1:idx2)))  % If No Ice Bottom, skip
            continue;
        end
        if isnan(mean(A.bt_elev_new(idx1:idx2)))  % If No Ice Bottom, skip
            continue;
        end
        
        p = polyfit(distance(idx1:idx2),A.bt_elev_new(idx1:idx2),1);    % Polygonal fitting
        slopeval(idx1:idx2) = polyval(p,distance(idx1:idx2));           % Ploygonal values after fitting
        base = distance(idx2)-distance(idx1);
        perpendicular = slopeval(idx2)-slopeval(idx1);
        hypotenuse = sqrt(base^2+perpendicular^2);
        angle(rline)=asin(perpendicular/hypotenuse)*180/pi;
        slopeerror = slopeval(idx1:idx2)-slopeval(idx1);         % Error betn original and fitting line
        dtime = 2*slopeerror/c/sqrt(er_ice);
        if rline==1
            figure(8);plot([idx1:idx2],A.bt_elev_new(idx1:idx2));
            hold on;plot([idx1:idx2],slopeval(idx1:idx2),'r--')
            figure(9);plot(10*log10(abs(A.Data_new(:,idx2)).^2))
        end
        for idx = idx1:idx2
            A.Data_new(:,idx) = interp1(elev_axis, A.Data_new(:,idx), elev_axis + slopeerror(idx-idx1 +1), 'linear',0);
            A.Elevation1_new(idx) = A.Elevation1_new(idx) - slopeerror(idx-idx1 +1);
            A.sf_elev_new(idx) = A.sf_elev_new(idx) - slopeerror(idx-idx1 +1);
            A.bt_elev_new(idx) = A.bt_elev_new(idx) - slopeerror(idx-idx1 +1);
            sf_new(idx) = sf_new(idx) + dtime(idx-idx1 +1);
            bt_new(idx) = bt_new(idx) + dtime(idx-idx1 +1);
        end
        if rline==1
            figure(8);hold on;plot([idx1:idx2],A.bt_elev_new(idx1:idx2),'g--');
            figure(9);hold on;plot(10*log10(abs(A.Data_new(:,idx2)).^2),'g--');
        end
    end
    
    
    fh = figure(3);
    %figure(fh);imagesc([1:size(A.Data_new,2)],elev_axis,10*log10(abs(A.Data_new).^2));title('Slope Correction Data')
    figure(fh);imagesc(distance,elev_axis,10*log10(abs(A.Data_new).^2));title('Slope Correction Data')
    ax = gca;
    ax.YDir = 'normal';
    hold on;plot(distance,A.sf_elev_new,'--');plot(distance,A.bt_elev_new,'--');
    
    
    
    %
     % Truncate data around ice bottom within bt.range_bins
      bt_range_bins =[-50:100];
      bt.val = NaN*ones(1,size(A.Data_new,2));
      bt.idx = NaN*ones(1,size(A.Data_new,2));
      bt.waveform = NaN*ones(length(bt_range_bins),size(A.Data_new,2));
      bt.inc_wf_ave=NaN*ones(size(bt.waveform,1),Nx);
      for rline = 1:size(A.Data_new,2)
        if ~isnan(A.bt_elev_new(rline)) & ~isinf(A.bt_elev_new(rline))
         % bt.idx(rline) = find(elev_axis<=A.bt_elev_new(rline),1,'first');
         bt.idx(rline) = round(interp1(elev_axis,[1:length(elev_axis)],A.bt_elev_new(rline)));
         bt.val(rline) = A.Data_new(bt.idx(rline),rline);
          first_idx = bt.idx(rline) + bt_range_bins(1);
          last_idx = bt.idx(rline) + bt_range_bins(end);
          if first_idx < 1 | last_idx>size(A.Data_new,1)
            bt.idx(rline) = NaN;
            bt.val(rline) = NaN;
            continue
          end
          lower=bt.idx(rline)+bt_range_bins(1);
          upper=bt.idx(rline)+bt_range_bins(end);
          if upper >size(A.Data_new,1)
              upper=size(A.Data_new,1);
          end
          bt.waveform(1:51+upper-bt.idx(rline),rline) = A.Data_new(lower:upper, rline);
        else
          continue
        end
      end
    
      debug_flag=0;
     if debug_flag==1
        figure(31);imagesc(lp(A.Data_new));
        hold on;plot(bt.idx,'--');
        figure(33);plot(lp(bt.waveform));
        figure(32);imagesc(lp(bt.waveform));
     end
    
    
    
    %}
    %% Coherene Index Calculation with Slope error Correction
    
    
%     for rline=1:size(A.Data_new,2)
%        if isnan(A.bt_elev_new(rline))
%            continue;
%        end
%         bt_idx_m = round(interp1(elev_axis,[1:length(elev_axis)],A.bt_elev_new));
%         [pwr,bt_idx]=max(power(bt_idx_m-5:bt_idx_m+5,rline));
%         A.bt_elev_new(rline)=(A.bt_elev_new(rline)-5*dr)+(bt_idx-1)*dr;
%     end
%     
%     A.bt_elev_new=sgolayfilt(A.bt_elev_new,1,Nx_int);
%     
%     fh=figure(100);
%     figure(fh);imagesc(distance,elev_axis,10*log10(abs(A.Data_new).^2));title('Slope Correction Data')
%     hold on;plot(distance,A.sf_elev_new,'--');plot(distance,A.bt_elev_new,'--');
%     ax = gca;
%     ax.YDir = 'normal';
    

    Nx0 = size(A.Dataweiner,2);
    Nx = floor(Nx0/Nx_int);
    Nx_mod = mod(Nx0,Nx_int);
    if Nx_mod>= Nx_int/2;
        Nx = Nx + 1;
    end
    square_int = zeros(size( A.Data_new,1),Nx); % incoherent integration, take square of abs first, then sum
    square_int_dB=zeros(size( A.Data_new,1),Nx);
    int_square = zeros(size(A.Data_new,1),Nx); % coherent integration, sum complex data first, then take square of abs
    coh_index_slopecorr =NaN*ones(1,Nx);
    Abruptiveindex=NaN*ones(1,Nx);
    Padj=NaN*ones(1,Nx);
    Latitude_mean=zeros(1,Nx);
    Longitude_mean=zeros(1,Nx);
    GPS_time_ave=zeros(1,Nx);
    
    risingedge=NaN*ones(1,Nx);
    fallingedge=NaN*ones(1,Nx);
    peakindex=NaN*ones(1,Nx);
    for rline =1:Nx
        idx1 = (rline-1)*Nx_int + 1;
        idx2 = rline*Nx_int;
        if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
            idx2 = Nx0;
        else
            idx2 = min(idx2,Nx0);
        end
        
        square_int(:,rline) = mean(abs(A.Data_new(:,idx1:idx2)).^2,2);
        int_square(:,rline) = abs(mean(A.Data_new(:,idx1:idx2),2)).^2;
        square_int_dB(:,rline) = 10*log10(square_int(:,rline));
        
        Latitude_mean(rline)=mean(A.Latitude(idx1:idx2));  %Mean Latitude
        Longitude_mean(rline)=mean(A.Longitude(idx1:idx2)); %Mean Longitude
        GPS_time_ave(rline)=mean(A.GPS_time(idx1:idx2));   %Mean GPS Time
        
        meanbt=nanmean(A.bt_elev_new(idx1:idx2));
        if meanbt==0 | isnan(meanbt)
            continue;               %skip if no ice bottom
        end
        
        
        meansf=nanmean(A.sf_elev_new(idx1:idx2));  %If bottom=surface; skip
        meanbt1=nanmean(bttime(idx1:idx2));
        if meanbt1==meansf
            continue;
        end
         depth=meansf-meanbt;
         if depth<150
             continue;
         end
        %bt_idx_m = find(elev_axis<=meanbt,1,'first');
        bt_idx_m = round(interp1(elev_axis,[1:length(elev_axis)],meanbt));
        
        if isempty(bt_idx_m)
            continue;  %looking for ice bottom
        end
        b1=bt_idx_m-5;               % Peak search index peak bins
        b2=bt_idx_m+5;
        if b2>size(A.Data_new,1)
            b2=size(A.Data_new,1);
        end
        
        [bt_val,bt_idx] = max(square_int(b1:b2,rline));  %Peak Index and Value
        bt_idx = bt_idx_m-5+bt_idx-1;
        bt_pwr = 10*log10(bt_val);
        
        if bt_pwr==0
            continue;
        end
        
        noise_bin1=bt_idx+270;
        noise_bin2=bt_idx+300;
        if noise_bin1>size(A.Data_new,1) || noise_bin2>size(A.Data_new,1)
            noise_bin1=size(A.Data_new,1)-30;
            noise_bin2=size(A.Data_new,1);
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
        if bt_idx-risingedge(rline)>50
            risingedge(rline)=bt_idx-50;  %Max 50 bins away
        end
        
        if bt_idx-risingedge(rline)<3
            risingedge(rline)=bt_idx-3;  %Guard band of 3
        end
        [r,tmp_idx]=(min(square_int_dB(risingedge(rline):bt_idx-3,rline)));
        risingedge=risingedge+tmp_idx-1;
        
        
        fallingedge(rline) = bt_idx+1;
        if fallingedge(rline)>size(A.Data_new,1)
            fallingedge(rline)=size(A.Data_new,1);
        end
        while square_int_dB(fallingedge(rline),rline)-noise > 0.05*SNR & fallingedge(rline)<size(A.Data_new,1) %Depth Bins risingedge(rline) and fallingedge(rline)
            fallingedge(rline) = fallingedge(rline) + 1;
        end
        if fallingedge(rline)-bt_idx>100
            fallingedge(rline)=bt_idx+100;
        end
        if fallingedge(rline)-bt_idx<3
            fallingedge(rline)=bt_idx+3;
        end
        coh_index_slopecorr(rline) = sum(int_square(risingedge(rline):fallingedge(rline),rline))/sum(square_int(risingedge(rline):fallingedge(rline),rline)); %Coherence Index
        
        Imeanx=sum((square_int(risingedge(rline):fallingedge(rline),rline)));  %Power summed from risingedge(rline) to fallingedge(rline)
        Abruptiveindex(rline)=bt_val/Imeanx;    %Abruptive Index
        peakindex(rline)=bt_idx;
      
        D1=risingedge(rline);
        D2=fallingedge(rline);
        Dnoise=noise;
        depth=meansf-meanbt;
        
        B=2.3*3000/(depth+2000);
        Padj(rline)=(lp(Imeanx*depth.^2)+B*depth/100)+15;
        
        bx1=bt_idx-50;
        bx2=bt_idx+100;
        if bx1<=0 |bx2>size(square_int_dB,1)
            continue;
        end
        bt.inc_wf_ave(:,rline) = square_int_dB(bt_idx-50:bt_idx+100,rline);
        
    end
    
%     Padj(Padj<(max(Padj)-40))=max(Padj)-40;
    distancenx=(max(distance).*Nx_int*0.5/(size(distance,2)):max(distance).*Nx_int/size(distance,2):max(distance).*Nx_int/size(distance,2)*Nx);
    hold on; figure(4); hold on; plot(distancenx,coh_index_slopecorr,'g','Displayname','SlopeCorrected'); title('Comparision of Coherence Index')
    legend('show');
    distancenx=distancenx/1000;
    figure(7);plot(distancenx,coh_index_slopecorr,'r','Displayname','Coherence Index'); title('Basal Roughness')
    figure(7); hold on; plot(distancenx,Abruptiveindex,'b','Displayname','Abruptive Index');
    xlabel('Distance(km)');
    legend('show');
    figure(6);plot(distancenx,Padj),title('Adjusted Intesnity')
    figure(11);plot(distancenx,coh_index_slopecorr,'r');ru
    figure(12); plot(distancenx,Abruptiveindex,'b');
    toc
    figure(4);plot(distancenx,coh_index_ogdata,'Displayname','Original Data');
   
    figure;subplot(3,1,1);plot(distancenx,Padj)
    subplot(3,1,2);plot(distancenx,coh_index_slopecorr)
    subplot(3,1,3);plot(distancenx,Abruptiveindex)
end