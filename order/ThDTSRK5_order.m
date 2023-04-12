function y_error = ThDTSRK5_order(a, dt, t, y0, nt)
f = @(t,y) a*(y-sin(t))+cos(t) ;
g = @(t,y) a*(a*(y-sin(t)) )-sin(t) ;
h= @(t,y) a*(a*(a*(y-sin(t)) ) ) - cos(t) ;

y = zeros(1,nt+1);
y(1) = y0;

%  ThDTSRK5_coefficients
tht=0.0;
cc1 = 0.198389107020261;
w1 = 0.501187671012343;
w2 = 0.167743974813318;
vv2 = 0.657963316199564;
ww2 = 1.49119408435601;
vvv2 = 0.119984650586874;
www2 = 0.0621952996182998;

cc2=cc1*cc1/2;
cc3=cc1*cc1*cc1/6;
a21=cc1;

vvv1=(1/60).*(23+(-60).*a21.*vv2+(-120).*a21.^2.*vv2+(-60).*a21.^3.* ...
    vv2+(-60).*vvv2+(-240).*a21.*vvv2+(-180).*a21.^2.*vvv2+(-5).*w1+( ...
    -5).*w2+30.*a21.^2.*w2+60.*a21.^3.*w2+60.*a21.^2.*ww2+(-60).* ...
    a21.^3.*ww2+120.*a21.*www2+(-180).*a21.^2.*www2);

vv1=(1/20).*(3+(-20).*vv2+60.*a21.^2.*vv2+40.*a21.^3.*vv2+120.*a21.* ...
    vvv2+120.*a21.^2.*vvv2+10.*w1+10.*w2+20.*a21.*w2+(-40).*a21.^3.* ...
    w2+(-60).*a21.^2.*ww2+40.*a21.^3.*ww2+(-120).*a21.*www2+120.* ...
    a21.^2.*www2);

ww1=(1/20).*(7+(-60).*a21.^2.*vv2+(-40).*a21.^3.*vv2+(-120).*a21.* ...
    vvv2+(-120).*a21.^2.*vvv2+10.*w1+10.*w2+(-20).*a21.*w2+40.* ...
    a21.^3.*w2+(-20).*ww2+60.*a21.^2.*ww2+(-40).*a21.^3.*ww2+120.* ...
    a21.*www2+(-120).*a21.^2.*www2);

www1=(1/60).*(8+(-60).*a21.^2.*vv2+(-60).*a21.^3.*vv2+(-120).*a21.* ...
    vvv2+(-180).*a21.^2.*vvv2+5.*w1+5.*w2+(-30).*a21.^2.*w2+60.* ...
    a21.^3.*w2+(-60).*a21.*ww2+120.*a21.^2.*ww2+(-60).*a21.^3.*ww2+( ...
    -60).*www2+240.*a21.*www2+(-180).*a21.^2.*www2);

v1 = 1 - w1;
v2 = -w2;

%  ThDRK5_coefficients
bbb1=1/16;
bbb2=5/48  ;
a21=2/5;
a212=a21*a21;  a213=a212*a21;

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
        y(i+1) = y(i) + dt*Lt0+ dt2/2* Ltt0 +  dt3*(bbb1*Lttt0 +bbb2*Lttt1);
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




