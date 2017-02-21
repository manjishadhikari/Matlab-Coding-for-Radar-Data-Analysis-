
a=A.Time;
b=A.Bottom;
[row,col]=size(A.Data);
Maxpower=zeros(1,col);
Noise=zeros(1,col);
Peaksnr=zeros(1,col);
fivepercentpower=zeros(1,col);
Imeanx=zeros(1,col);
Abruptiveindex=zeros(1,col);

d=zeros(1,200);
ko=zeros(1,200);

for i=1:col
   
tmp1=abs(a-b(i));
[r,c]=min(tmp1);
Maxpower(i)=A.Data(c,i);
Noise(i)=mean(A.Data(c+200:1600,i));
Peaksnr(i)=Maxpower(i)-Noise(i);
fivepercentpower(i)=Noise(i)+0.05*Peaksnr(i);
tmp2=(abs(A.Data(:,i)-fivepercentpower(i)));
tmp3=sort(tmp2);

    for j=1:200
       ko=find(tmp2==tmp3(j));   
       d(1,j)=ko(1);
    end
        G=d(c>d);
            tmp4=abs(c-G);
            tmp5=sort(tmp4);
            val1=tmp5(1);
            l=find(tmp4==val1);
        H=d(c<d);
            tmp6=abs(c-H);
            tmp7=sort(tmp6);
            val2=tmp7(1);
            u=find(tmp6==val2);  
h1=G(l(1));
h2=H(u(1));
        if h1<h2
            Imeanx(i)=sum(A.Data(h1:h2,i));
        elseif h1>h2
            Imeanx(i)=sum(A.Data(h2:h1,i));
        end

 Abruptiveindex(i)=Maxpower(i)/Imeanx(i); 
end
figure
plot(Abruptiveindex)
title('Abruptive Index Profile');

ylabel('Abruptive Index');
    

