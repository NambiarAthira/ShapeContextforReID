
function [mini]=filterr(vari)
% y =[23.7111   21.3362   19.1567   14.9394   13.4477   13.7501   12.6787   13.0537   16.4223   17.6788   12.4185    9.5504    5.3949    8.2945 5.8001    6.8703   11.8432   11.7708   14.7936   19.7456    7.7598    5.0385    7.2796    7.9664   11.6833   13.8809   16.9107   17.0069 16.6278   14.3378    9.5948    8.3242    4.4811    4.1053    4.4120    6.4427]
y=vari;
% y =[20.5976   13.5701   13.9627    9.9663    7.9642    6.7409   14.9115   12.6017   10.0093    9.0040    7.5536    8.8529    7.4730 6.4695    6.1232   10.9862   11.8515   14.7389   14.5416   13.3809   13.2495   11.4038   11.8744    8.6399    4.8503    6.5282 5.9251    6.7610    6.1624    6.3913    6.4917]
x=1:1:size(y,2);
% figure,plot(x,y,'b')
% hold on
% a = 1;
% b = [1/6 1/6 1/6 1/6];
% y = filter(b,a,x)

% figure
% subplot(2,1,1)
% plot(x,y)
[ymax,imax,ymin,imin] = extrema(y);
% hold on
% plot(x(imax),ymax,'k*',x(imin),ymin,'g*')
output = tsmovavg(y, 's', 3);
output=smooth(y,'moving');
% hold on,plot(output,'r--')
% xlabel('Row number of the image');
% ylabel('pixel counts in each row');
% legend('raw data','maxima','minima','filtered data')
% title('Initial measurements')
% subplot(2,1,2)
% plot(x,output)
% xlabel('Row number of the image');
% ylabel('pixel counts in each row');
% title('Filtered output of the point Moving Average filter')
[ymax2,imax2,ymin2,imin2] = extrema(output);
% hold on
% plot(x(imax2),ymax2,'m*',x(imin2),ymin2,'g*')
% legend('filtered data','maxima','minima')
mini=imin2;
if size(mini)==1
    mini=round(vertcat(mini,(imin+imax)/2));
end
mini;

% subplot(3,1,3)
% peak=zeros(size(output));
% max=sort(x(imax2))
% x=find(diff(max)<5) 
% max(x)=(max(x)+max(x+1))/2
% max(x+1)=[]
% max=round(max)
% peak(max)=1
% plot(peak)
% title('The positions of peaks')
% xlabel('Frame Number');
% ylabel('Extrema');
% saveas(gcf,'gaitcycle.fig')
end

% fid = fopen('vari.txt', 'wt');
% fprintf(fid, [repmat('%g\t', 1, size(vari,2)-1) '%g\n'], vari.');
% fclose(fid)


% tdata = (1:length(vari))'; 
% p_coeffs = polyfit(tdata,vari,6);
% 
% figure
% plot(vari,'o-')
% hold on
% tfit = (1:0.01:length(vari))';
% yfit = [ones(size(tfit)) cos((2*pi/12)*(tfit-7))]*s_coeffs;
% plot(tfit,yfit,'r-','LineWidth',2)