
%% Determination of coherence Index 

fc = 195e6;
fs = 1.1111e8;  %Sampling Frequency of the radar
c = 3e8;        %Speed of light
er_ice = 3.15;  %Permittivity of ice
load_flag = 1;
if load_flag
     A = load('/cresis/snfs1/dataproducts/ct_data/rds/2011_Antarctica_TO/CSARP_ant/20111201_04/Data_20111201_04_002.mat');
    L = load('/cresis/snfs1/dataproducts/ct_data/rds/2011_Antarctica_TO/CSARP_post/CSARP_layerData/20111201_04/Data_20111201_04_002.mat');
end
sf = interp1(L.GPS_time,L.layerData{1}.value{2}.data,A.GPS_time); %Two way prop time to surface
bt = interp1(L.GPS_time,L.layerData{2}.value{2}.data,A.GPS_time); %Two way prop time to bottom
figure(10); imagesc([],A.Time,10*log10(abs(A.Data)));
hold on;plot(sf,'--');plot(bt,'--');

%% Correction for Elevation
sf_new = zeros(size(sf));
bt_new = zeros(size(bt));
A.Elevation_new = zeros(size(A.Elevation));
A.sf_elev_new = zeros(size(A.Elevation));
A.bt_elev_new = zeros(size(A.Elevation));

% 1)Remove data before zero time
negative_bins = A.Time < 0;
A.Time_new = A.Time(~negative_bins);
A.Data_new = A.Data(~negative_bins,:);

% 2)Create elevation axis to interpolate to
[max_elev,max_elev_idx] = max(A.Elevation);
min_elev = min(A.Elevation - sf*c/2 - (A.Time_new(end)-sf)*c/2/sqrt(er_ice));
dt = A.Time(2)-A.Time(1);
dr = dt * c/2 / sqrt(er_ice);
dt_air = dr/(c/2);
elev_axis = max_elev:-dr:min_elev;
new_time = zeros(length(elev_axis),length(A.Elevation));

% 3)Zero pad data to create space for interpolated data
zero_pad_len = length(elev_axis) - length(A.Time_new);
A.Data_new = cat(1,A.Data_new,zeros(zero_pad_len,size(A.Data_new,2)));

% 4)Determine the corrections to apply to elevation and layers
dRange = max_elev - A.Elevation;
dBins = round(dRange / (c/2) / dt);
dtime = dRange/(c/2);

for rline = 1:size(A.Data_new,2)
    % Determine elevation bins before surface
    sf_elev = A.Elevation(rline) - sf(rline) * c/2;
    time0 = -(max_elev - A.Elevation(rline))/(c/2);
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
    A.Elevation_new(rline) = A.Elevation(rline) + dRange(rline);
    sf_new(rline) = sf(rline) + dtime(rline);
    bt_new(rline) = bt(rline) + dtime(rline);
    A.sf_elev_new(rline) = A.Elevation_new(rline) - sf_new(rline)*c/2;
    A.bt_elev_new(rline) = A.sf_elev_new(rline) - (bt_new(rline)-sf_new(rline))*c/2/sqrt(er_ice);
end
fh = figure(1);
figure(fh);imagesc([],elev_axis,10*log10(abs(A.Data_new).^2));title('Elevation Correction Data')
ax = gca;
ax.YDir = 'normal';
hold on;plot(A.sf_elev_new,'--');plot(A.bt_elev_new,'--');
bt_slope = diff(A.bt_elev_new)./diff(distance);
figure(2);plot(10*log10(abs(A.Data_new(:,6081:6112)).^2)); title('Before Slope Correction');


%% Ice Slope Correction

Nx_int = 32;
Nx0 = size(ElevCorrectedData,2);
Nx = floor(Nx0/Nx_int);
Nx_mod = mod(Nx0,Nx_int);
if Nx_mod>= Nx_int/2;
    Nx = Nx + 1;
end
slopeerror=zeros(1,Nx);
slopeval=zeros(1,Nx);
angle=zeros(1,Nx);
for rline = 1:Nx
    
    idx1 = (rline-1)*Nx_int + 1;
    idx2 = rline*Nx_int;
    if Nx0 - idx2 > 0 && Nx0 - idx2 < Nx_int/2;
        idx2 = Nx0;
    else
        idx2 = min(idx2,Nx0);
    end
    
    p = polyfit(distance(idx1:idx2),A.bt_elev_new(idx1:idx2),1);    % Polygonal fitting
    slopeval(idx1:idx2) = polyval(p,distance(idx1:idx2));           % Ploygonal values after fitting
    base = distance(idx2)-distance(idx1);
    perpendicular = slopeval(idx2)-slopeval(idx1);
    hypotenuse = sqrt(base^2+perpendicular^2);
    angle(rline)=asin(perpendicular/hypotenuse)*180/pi;
    slopeerror = slopeval(idx1:idx2)-slopeval(idx1);         % Error betn original and fitting line
    dtime = 2*slopeerror/c/sqrt(er_ice);
    if idx1 == 6081
        figure(5);plot([idx1:idx2],A.bt_elev_new(idx1:idx2));
        hold on;plot([idx1:idx2],slopeval(idx1:idx2),'r--')
        figure(6);plot(10*log10(abs(A.Data_new(:,idx2)).^2))
    end
    for idx = idx1:idx2
        A.Data_new(:,idx) = interp1(elev_axis, A.Data_new(:,idx), elev_axis + slopeerror(idx-idx1 +1), 'linear',0);
        A.Elevation_new(idx) = A.Elevation_new(idx) - slopeerror(idx-idx1 +1);
        A.sf_elev_new(idx) = A.sf_elev_new(idx) - slopeerror(idx-idx1 +1);
        A.bt_elev_new(idx) = A.bt_elev_new(idx) - slopeerror(idx-idx1 +1);  
        sf_new(idx) = sf_new(idx) + dtime(idx-idx1 +1);
        bt_new(idx) = bt_new(idx) + dtime(idx-idx1 +1);
    end
    if idx1 == 6081
        figure(5);hold on;plot([idx1:idx2],A.bt_elev_new(idx1:idx2),'g--');
        figure(6);hold on;plot(10*log10(abs(A.Data_new(:,idx2)).^2),'g--');
    end
end
fh = figure(3);
figure(fh);imagesc([],elev_axis,10*log10(abs(A.Data_new).^2));title('Slope Correction Data')
ax = gca;
ax.YDir = 'normal';
hold on;plot(A.sf_elev_new,'--');plot(A.bt_elev_new,'--');
figure(4);plot(10*log10(abs(A.Data_new(:,6081:6112)).^2)); title('After Slope Correction')






















