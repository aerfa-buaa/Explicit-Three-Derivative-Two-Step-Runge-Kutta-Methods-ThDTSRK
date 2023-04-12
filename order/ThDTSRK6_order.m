function y_error = ThDTSRK6_order(a, dt, t, y0, nt)
f = @(t,y) a*(y-sin(t))+cos(t) ;
g = @(t,y) a*(a*(y-sin(t)) )-sin(t) ;
h= @(t,y) a*(a*(a*(y-sin(t)) ) ) - cos(t) ;

y = zeros(1,nt+1);
y(1) = y0;

%  ThDTSRK6_coefficients
  tht=0.0;
cc1 = 0.587325896573798;
 ww2 = -0.141169152307059;
 vvv2 = 0.0243580486114999;
www2 = -0.0227607642077618;

cc2=cc1*cc1/2;
cc3=cc1*cc1*cc1/6;
a21=cc1;

vvv1 =(1/120).*(111+tht+(-120).*vvv2+(-1080).*a21.*vvv2+(-2160).* ...
  a21.^2.*vvv2+(-1200).*a21.^3.*vvv2+120.*a21.*ww2+360.*a21.^2.*ww2+ ...
  1200.*a21.^3.*ww2+(-360).*a21.*www2+1440.*a21.^2.*www2+(-1200).* a21.^3.*www2);

vv1 =(1/10).*((-31)+(-1).*tht+360.*a21.*vvv2+960.*a21.^2.*vvv2+600.* ...
  a21.^3.*vvv2+10.*ww2+(-60).*a21.^2.*ww2+(-600).*a21.^3.*ww2+240.* ...
  a21.*www2+(-840).*a21.^2.*www2+600.*a21.^3.*www2);

v1 =(1/2).*(15+tht+(-120).*a21.*vvv2+(-360).*a21.^2.*vvv2+(-240).* ...
  a21.^3.*vvv2+240.*a21.^3.*ww2+(-120).*a21.*www2+360.*a21.^2.*www2+ ...
  (-240).*a21.^3.*www2);
w1=(1/2).*((-13)+tht+120.*a21.*vvv2+360.*a21.^2.*vvv2+240.*a21.^3.* ...
  vvv2+(-240).*a21.^3.*ww2+120.*a21.*www2+(-360).*a21.^2.*www2+240.* ...
  a21.^3.*www2);

ww1=(1/10).*((-29)+tht+240.*a21.*vvv2+840.*a21.^2.*vvv2+600.*a21.^3.* ...
  vvv2+(-10).*ww2+60.*a21.^2.*ww2+(-600).*a21.^3.*ww2+360.*a21.* ...
  www2+(-960).*a21.^2.*www2+600.*a21.^3.*www2);

www1= (1/120).*((-49)+tht+360.*a21.*vvv2+1440.*a21.^2.*vvv2+1200.* ...
  a21.^3.*vvv2+(-120).*a21.*ww2+360.*a21.^2.*ww2+(-1200).*a21.^3.* ...
  ww2+(-120).*www2+1080.*a21.*www2+(-2160).*a21.^2.*www2+1200.*  a21.^3.*www2);
 v2 = 0;  w2 = 0;  vv2 = -ww2;

%  ThDRK7_coefficients
bbb1=1/30;
bbb2=1/15  +13*sqrt(2)/480  ;
bbb3=1/15  - 13*sqrt(2)/480  ;
a21=3/7- sqrt(2)/7;
a31=3/7+ sqrt(2)/7;
aaa32=122/7203+ 71*sqrt(2)/7203;

a212=a21*a21;  a213=a212*a21;
a312=a31*a31;  a313=a312*a31;
for i=1:nt+1
    if  i<3
        dt=dt/2;
        dt2=dt*dt; dt3=dt2*dt;
        Lt0 = f((i-1)*dt,y(i));  Ltt0 = g((i-1)*dt,y(i));   Lttt0 = h((i-1)*dt,y(i));
        y1=y(i)+a21* dt*Lt0 + a212 /2* dt2*Ltt0 +a213 /6* dt3*Lttt0 ;
        t1 = (i-1)*dt +a21* dt;
        Lt1  = f( t1 ,  y1 );
        Ltt1 = g( t1, y1 );
        Lttt1 = h( t1, y1 );
        y2=y(i)+a31* dt*Lt0+ a312 /2* dt2*Ltt0 + (a313 /6 - aaa32) * dt3*Lttt0 + aaa32 * dt3*Lttt1;
        t2 = (i-1)*dt +a31* dt;
        Lt2  = f( t2,  y2  );
        Ltt2  = g( t2 , y2   );
        Lttt2  = h( t2, y2  );
        y(i+1) = y(i) + dt*Lt0+ dt2/2* Ltt0 +  dt3*(bbb1*Lttt0 +bbb2*Lttt1+bbb3*Lttt2);
        dt=dt*2;
        if i==1
            ii=1;
            dt2=dt*dt; dt3=dt2*dt;
            Lt0 = f((ii-1)*dt,y(i));  Ltt0 = g((ii-1)*dt,y(i));   Lttt0 = h((ii-1)*dt,y(i));
            y1=y(i)+cc1* dt*Lt0 + cc2* dt2*Ltt0 +cc3* dt3*Lttt0 ;
            t1 = (ii-1)*dt +cc1* dt;
            Lt1  = f( t1 ,  y1 );            Ltt1 = g( t1, y1 );            Lttt1 = h( t1, y1 );
            yf=y(i);   Ltf0=Lt0;   Ltf1=Lt1;  Lttf0=Ltt0;  Lttf1=Ltt1; Ltttf0=Lttt0;  Ltttf1=Lttt1;
        end
    else
        ii=i-1;
        dt2=dt*dt; dt3=dt2*dt;
        Lt0 = f((ii-1)*dt,y(i));  Ltt0 = g((ii-1)*dt,y(i));   Lttt0 = h((ii-1)*dt,y(i));
        y1=y(i)+cc1* dt*Lt0 + cc2* dt2*Ltt0 +cc3* dt3*Lttt0 ;
        t1 = (ii-1)*dt +cc1* dt;
        Lt1  = f( t1 ,  y1 );
        Ltt1 = g( t1, y1 );
        Lttt1 = h( t1, y1 );
        y(i+1) =(1-tht) *y(i) +tht*yf + dt*(v1*Lt0 +v2*Lt1) + dt2*(vv1*Ltt0 +vv2*Ltt1) + dt3*(vvv1*Lttt0 +vvv2*Lttt1) ...
            + dt*(w1*Ltf0 +w2*Ltf1) + dt2*(ww1*Lttf0 +ww2*Lttf1) + dt3*(www1*Ltttf0 +www2*Ltttf1);
        yf=y(i);   Ltf0=Lt0;   Ltf1=Lt1;  Lttf0=Ltt0;  Lttf1=Ltt1; Ltttf0=Lttt0;  Ltttf1=Lttt1;
    end
end
y_error=abs (sin(t)-y(i+1) );
end




