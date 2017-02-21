close all
load_flag = 1;
if load_flag
     A = load('/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_post/CSARP_qlook/20120327_01/Data_20120327_01_001.mat');
    L = load('/cresis/snfs1/dataproducts/public/data/rds/2012_Greenland_P3/CSARP_layerData/20120327_01/Data_20120327_01_001.mat');
end
dTime = A.Time(2)-A.Time(1);
bt = interp1(L.GPS_time,L.layerData{2}.value{2}.data,A.GPS_time);
Nx_int = 32;
Nx0 = size(A.Data,2);
Nx = floor(Nx0/Nx_int);
Nx_mod = mod(Nx0,Nx_int);
if Nx_mod>= Nx_int/2;
    Nx = Nx + 1;
end
square_int = zeros(size(A.Data,1),Nx); % incoherent integration, take square of abs first, then sum 
int_square = zeros(size(A.Data,1),Nx); % coherent integration, sum complex data first, then take square of abs
coh_index = zeros(1,Nx);

for rline = 1:Nx
    idx1 = (rline-1)*Nx_int + 1;
    idx2 = rline*Nx_int;
    if Nx0 - idx2 > 0 & Nx0 - idx2 < Nx_int/2;
        idx2 = Nx0;
    else
        idx2 = min(idx2,Nx0);
    end
    square_int(:,rline) = mean(abs(A.Data(:,idx1:idx2)).^2,2);
    int_square(:,rline) = abs(mean(A.Data(:,idx1:idx2),2)).^2;   
    bt_idx_m = find(A.Time>mean(bt(idx1:idx2)),1,'first');
    [bt_val,bt_idx] = max(square_int(bt_idx_m-50:bt_idx_m+50,rline));
    bt_idx = bt_idx + bt_idx_m -50 - 1; 
    bt_pwr = 10*log10(bt_val);
    noise_bin1 = bt_idx+500;
    noise_bin2 = bt_idx+530;
    noise = 10*log10(mean(mean(abs(A.Data(noise_bin1:noise_bin2,:)).^2)));
    D1 = bt_idx-1;
    SNR = bt_pwr-noise;
    square_int_dB = 10*log10(square_int(:,rline));
    % find D1 and D2 for integration in range
    while square_int_dB(D1)-noise > 0.05*SNR & bt_idx-D1<50
        D1 = D1 - 1;
    end
    D2 = bt_idx+1;
    while square_int_dB(D2)-noise > 0.05*SNR
        D2 = D2 + 1;
    end
    coh_index(rline) = sum(int_square(D1:D2,rline))/sum(square_int(D1:D2,rline));
end
figure(1);plot(coh_index);
