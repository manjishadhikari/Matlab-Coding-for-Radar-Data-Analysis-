# Matlab-Coding-for-Radar-Data-Analysis-
%These are the radar data processing projects that I have been currently working on.
%Loadingdata

A=load('Data_20110506_01_027.mat');
B=load('Layers_20110506_01_027.mat');

subplot(2,1,1)
imagesc(A.GPS_time,A.Time,10*log10(A.Data)) %Imageplot of Icebed "Q.no 1"
xlabel('GPS time (sec)');
hc = colorbar;
set(get(hc,'ylabel'),'string','Relative power (dB)');
title('Radar iamging of ice sheet');
hold on;


%Icepicklayers
for i=1:346
y=[0 B.layerData{i}.value{2}.data];
plot(A.GPS_time,y,'Color','w');  %Plottingicepicklayers "Q.No. 2"
end


[row, column]=size(A.GPS_time);
C=[A.GPS_time zeros(1,3000-column)];
F=reshape(C,[100,30]);
P=mean(F); %Mean of every 100 GPS_time points "Q. NO 3"


for j=1:346
    hold on
d=[0 B.layerData{j}.value{2}.data];
for i=1:30
in=interp1(A.GPS_time,d,P); %Two way prop time by interpolation   "Q.No 3"
plot(P,in,'b*')  %Plotting the  ice layers at mean GPS points  "Q. No.4"
hold on
end
end
