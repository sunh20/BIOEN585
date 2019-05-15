clear all; close all; 
% Data are from 
% http://www-rohan.sdsu.edu/~jmahaffy/courses/f00/math122/labs/labj/lotvol.xls
% Year   Hares (x1000)   Lynx(x1000)
global data
data=...
    [1900   30  4
    1901    47.2    6.1
    1902    70.2    9.8
    1903    77.4    35.2
    1904    36.3    59.4
    1905    20.6    41.7
    1906    18.1    19
    1907    21.4    13
    1908    22  8.3
    1909    25.4    9.1
    1910    27.1    7.4
    1911    40.3    8
    1912    57  12.3
    1913    76.6    19.5
    1914    52.3    45.7
    1915    19.5    51.1
    1916    11.2    29.7
    1917    7.6 15.8
    1918    14.6    9.7
    1919    16.2    10.1
    1920    24.7    8.6];
% guess the parameters [r a b m] and initial conditions[hares lynx]:
guesses = [1.00 0.02 0.02 1.00 35.7 3.99];
% optimize the parameters and initial conditions:
[estimates, J] = fminsearch(@obj,guesses);
disp(estimates);
% graph the model behavior using these optimized values:
tspan = data(:,1)-1900;
y0 = [estimates(5) estimates(6)];
options=[];
[t,y] = ode23s(@lotkavolterraODE,tspan,y0,options,estimates);
plot(t,data(:,2),'o',t,y(:,1),t,data(:,3),'*',t,y(:,2))



