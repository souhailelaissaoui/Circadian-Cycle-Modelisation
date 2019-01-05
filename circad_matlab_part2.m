
% 14/11/2017

%% Parameters
clear all
close all

v1=0.7 ;
v2 = 0.35 ;
v4 = v2 ;
v6 = v2 ;

K1 = 1 ;
K2 = K1 ;
K4 = K1 ;
K6 = K1 ;

k3 = 0.7 ;
k5 =k3 ;


%% SIMULATION OF THE MODEL

%% Question 4
%trajectoires  - Euler
dt1=0.01;
n1=10000/dt1;
x=zeros(1,n1);
y=zeros(1,n1);
z=zeros(1,n1);
x(1)=0; %x0
y(1)=1; %y0
z(1)=5; %z0

for i=1:n1-1
   x(i+1)=x(i)+dt1*(v1 * K1^4/(K1^4+z(i)^4)- v2 * x(i)/(K2+x(i))) ;  
   y(i+1)=y(i)+dt1*(k3 * x(i)- v4 * y(i)/(K4+y(i))); 
   z(i+1)=z(i)+dt1*(k5 * y(i)- v6 * z(i)/(K6+z(i))); 
end

figure
subplot(3,1,1)
plot(x)

subplot(3,1,2)
plot(y)

subplot(3,1,3)
plot(z)

figure
plot3(x,y,z); hold on 


%% Question 5

figure
[ppx, f] =periodogram(x(50000:length(x)),[],[],1/dt1) ;
plot(1./f, ppx)

figure
[ppy, f] =periodogram(y(50000:length(y)),[],[],1/dt1) ;
plot(1./f, ppy)

figure
[ppz, f] =periodogram(z(50000:length(z)),[],[],1/dt1) ;
plot(1./f, ppz)
%xlim([0 200])

%% Question 6

v1=linspace(0.001, 1, 50) ;
x_min=zeros(1,length(v1));
x_max=zeros(1,length(v1));
figure

for j=1:length(v1)
    x=zeros(1,n1);
    y=zeros(1,n1);
    z=zeros(1,n1);
    x(1)=0; %x0
    y(1)=1; %y0
    z(1)=5; %z0

    for i=1:n1-1
       x(i+1)=x(i)+dt1*(v1(j) * K1^4/(K1^4+z(i)^4)- v2 * x(i)/(K2+x(i))) ;  
       y(i+1)=y(i)+dt1*(k3 * x(i)- v4 * y(i)/(K4+y(i))); 
       z(i+1)=z(i)+dt1*(k5 * y(i)- v6 * z(i)/(K6+z(i))); 
    end

    plot3(x,y,z); hold on 
    
    x_min(j)=min(x) ;
    x_max(j)=max(x) ;
end

%% Question 7
figure
plot(v1,x_min, v1, x_max)

%% Question 8
v1=linspace(0.3, 1.3, 20) ;
period=zeros(1,length(v1));

for j=1:length(v1)
    x=zeros(1,n1);
    y=zeros(1,n1);
    z=zeros(1,n1);
    x(1)=0; %x0
    y(1)=1; %y0
    z(1)=5; %z0

    for i=1:n1-1
       x(i+1)=x(i)+dt1*(v1(j) * K1^4/(K1^4+z(i)^4)- v2 * x(i)/(K2+x(i))) ;  
       y(i+1)=y(i)+dt1*(k3 * x(i)- v4 * y(i)/(K4+y(i))); 
       z(i+1)=z(i)+dt1*(k5 * y(i)- v6 * z(i)/(K6+z(i))); 
    end
    [ppx, f] = periodogram(x(50000:length(x))-mean(x),[],[],1/dt1) ;
    [maxi,index] = max(ppx) ;
    period(j) = 1./f(index) ;
end

plot(v1, period)

