

% function [cp,lp]= get_points(x)
I = imread('SC4_2.png')
BW2 = edge(I,'prewitt')
imshow(BW2)
[x,y]=find(BW2==1)
N=length(x)
% aux contem userdata de stopbutton

%imshow(x); 
% hold on,
% N=0; but=1;

% w = waitforbuttonpress;
% if w == 0
%     disp('Button press')
% else
%     disp('Key press')
% end
cp=[]; lp=[];

% % while (but ==1 | but ==2 | but ==3 | but == 32)
%     %pode acontecer q o médico queira acabar o contorno com "stop"
%     %     aux = get(handles.stopButton,'userdata');
%     %     if aux
%     %         set(handles.obsText,'string','Deve terminar contorno com ENTER');
%     %         return;
%     %     end
% 
%     if N<0, return, end
%     %plot(real(z_cR),imag(z_cR),'r-'); hold on,
%     [ci,li,but]=ginput(1);
%     %aux = get(handles.stopButton,'userdata');
%     %if aux
%     %    cp=[];lp=[];N=0;
%     %    return;
%     %end
%     %but
%     if but == 1 %add point

for i=1:N
        cp =  x(i);
        lp=  y(i);
        plot(cp,lp,'r.','MarkerSize',8); drawnow;hold on;
        if N > 1
            plot(cp,lp,'r.-','MarkerSize',8); drawnow;hold on;
        end
end
%     if but == 32 %remove point
%         hold off,
%         imshow(x); hold on,
%         plot(real(b_f)*1 ,imag(b_f)*1,'y.');  %drawnow
%         cp = cp(1:N-1);
%         lp = lp(1:N-1);
%         N  = N-1;
%         %N
%         plot(cp(:),lp(:),'r.-','MarkerSize',8); drawnow;
%     end
% end
cp=cp'; lp=lp';
% if but == []
%     return;
% end
% hold off;


