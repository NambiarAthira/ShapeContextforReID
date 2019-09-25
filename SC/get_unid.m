

function [cp,lp,NP]= get_unid(x);

N=0;

but=1;
while but ~= 3
    [ci,li,but]=ginput(1);  %add point 
    if but == 1             
        N=N+1;
        cp(N,1)=ci;
        lp(N,1)=li;
        plot(cp,lp,'b.','MarkerSize',8); drawnow;
        but;
    end
    if but == 2             %remove point
        %imshow(x); hold on,
        cp = cp(1:N-1);
        lp = lp(1:N-1);
        plot(cp(:),lp(:),'b.','MarkerSize',8); drawnow;
        %hold off,
        N = length(cp)
        but; 
    end
end
NP=N;

%hold off;  %descomentar depois

      
