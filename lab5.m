clc;clear;clf;
increment=.01;
cw=2;
dw=2.01;
SmoothWidth=4;

x=0:increment:20;
y=zeros(size(x));
y(900:1300)=1.3;

y=y+.01.*randn(size(y));

c=gaussian(x,0,cw)+gaussian(x,max(x),cw);
c2=gaussian(x,0,dw)+gaussian(x,max(x),dw);

yc=ifft(fft(y).*fft(c))./sum(c);

yc=yc+.00000001.*randn(size(yc));

ydc=ifft(fft(yc)./fft(c2)).*sum(c2);

subplot(2,2,1); plot(x,y); title('original y');
subplot(2,2,2); plot(x,yc(1:2001)); title(['yc=y convoluted with c. Width = ' num2str(cw) ]);
subplot(2,2,3); plot(x,ydc);title(['ydc=recovered y. Width = ' num2str(dw) ]);
subplot(2,2,4); plot(x,fastsmooth(ydc,SmoothWidth,3));title('smoothed recovered y');

function g = gaussian(x, mu, sigma)
    g = exp(-((x - mu).^2) / (2 * sigma^2));
end

function s = fastsmooth(a,w,type)
    s = movmean(a,w);
end