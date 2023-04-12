%%  ThDTSRK7_order
clc
clear all
t=double(2.8*pi);   % 
a=double(-10) ;
y0=0.0 ;
n_dt=[90, 125, 180, 250, 355, 500];
for i=1:6
    n_dt(i)
    dx(i)=t/n_dt(i);
    dt=dx(i);
    nt=round(t/ dt);
    tt=nt* dt
    error(i)=  ThDTSRK7_order(a, dt, tt, y0, nt) ;
    if i>1
        e_order1(i)=-log10(error(i)/error(i-1))/log10(dx(i)/dx(i-1));
    end
    errlog(i)= log10( error(i) );
    e_x(i)=log10(dx(i));
end


%%  ThDTSRK6_order
clc
clear all

t=double(2.8*pi);   % 
a=double(-10) ;
y0=0.0 ;
n_dt=[ 125, 180, 250, 355, 500, 700 ];
for i=1:6
    n_dt(i)
    dx(i)=t/n_dt(i);
    dt=dx(i);
    nt=round(t/ dt);
    tt=nt* dt
    error(i)=  ThDTSRK6_order(a, dt, tt, y0, nt) ;
    if i>1
        e_order1(i)=-log10(error(i)/error(i-1))/log10(dx(i)/dx(i-1));
    end
    errlog(i)= log10( error(i) );
    e_x(i)=log10(dx(i));
end

%%  ThDTSRK5_order
clc
clear all

t=double(2.8*pi);   % 
a=double(-10) ;
y0=0.0 ;
n_dt=[ 180, 250, 355, 500, 700, 1000];
for i=1:6
    n_dt(i)
    dx(i)=t/n_dt(i);
    dt=dx(i);
    nt=round(t/ dt);
    tt=nt* dt
    error(i)=  ThDTSRK5_order(a, dt, tt, y0, nt) ;
    if i>1
        e_order1(i)=-log10(error(i)/error(i-1))/log10(dx(i)/dx(i-1));
    end
    errlog(i)= log10( error(i) );
    e_x(i)=log10(dx(i));
end

 

