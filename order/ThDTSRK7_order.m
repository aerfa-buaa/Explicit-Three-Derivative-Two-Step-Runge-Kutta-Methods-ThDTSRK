function y_error = ThDTSRK7_order(a, dt, t, y0, nt)
f = @(t,y) a*(y-sin(t))+cos(t) ;
g = @(t,y) a*(a*(y-sin(t)) )-sin(t) ;
h= @(t,y) a*(a*(a*(y-sin(t)) ) ) - cos(t) ;

y = zeros(1,nt+1);
y(1) = y0;

%  ThDTSRK7_coefficients
cc1=0.5;  tht=0.0;
cc2=cc1*cc1/2;
cc3=cc1*cc1*cc1/6;
a21=cc1;
vvv1=(1/1680) *(a21+10 *a21 ^3) ^(-1) *((-209)+300 *a21+(-6270) *  ...
    a21 ^2+15540 *a21 ^3+tht+20 *a21 *tht+30 *a21 ^2 *tht+140 * a21 ^3 *tht);
vvv2=(1/1680) *(a21+10 *a21 ^3) ^(-1) *(209+(-1) *tht);
vv1=(1/28) *(1+10 *a21 ^2) ^(-1) *((-45)+627 *a21+(-868) *a21 ^2+(-3) *tht+(-3) *a21 *tht+(-28) *a21 ^2 *tht);
vv2=0;
v1=(1/14) *(1+10 *a21 ^2) ^(-1) *(105+(-627) *a21+1050 *a21 ^2+7 *tht+3 *a21 *tht+70 *a21 ^2 *tht);
v2=0 ;
www2=(1/1680) *(a21+10 *a21 ^3) ^(-1) *((-209)+tht);
www1=(1/1680) *(a21+10 *a21 ^3) ^(-1) *(209+(-1940) *a21+6270 *a21 ^2+...
    (-6860) *a21 ^3+(-1) *tht+20 *a21 *tht+(-30) *a21 ^2 *tht+140 * a21 ^3 *tht);
ww2=0 ;
ww1=(1/28) *(1+10 *a21 ^2) ^(-1) *((-123)+627 *a21+(-812) *a21 ^2+3 *tht+(-3) *a21 *tht+28 *a21 ^2 *tht);
w1=(1/14) *(1+10 *a21 ^2) ^(-1) *((-91)+627 *a21+(-910) *a21 ^2+7 * tht+(-3) *a21 *tht+70 *a21 ^2 *tht);
w2=0 ;

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




