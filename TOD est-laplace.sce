clc;
clear;
clf(1);

//generate theoritical data for small n
Ts=24;
Z=0.02;
p=0.6;
td=-2;
Td=38;
t=0:1.5:15
T=Ts+(Td-Ts)*((p/(p-Z))*exp(-Z*(t-td))-(Z/(p-Z)*exp(-p*(t-td))));
theta=T-Ts;
//disp(theta)

//find a and b by (3.6)
// ambil n=10
//disp(theta(10));
//disp(theta(9));
// compute b
b=log(theta($-1)/theta($))/(t($)-t($-1));
// compute a
a=theta($)*exp(b*t($));
printf("a = %f\n" ,a)
printf("b = %f\n" ,b)


t2=16.5:3:20000
t1=[t t2]
theta2=a*exp(-b*t2)
theta1=[theta theta2]
//compute the improper integral in (3.5)
r0 = intsplin(t1,theta1)
//disp(r0)
r1 = (-1)*intsplin(t1,theta1.*t1)
//disp(r1)
r2 = intsplin(t1,theta1.*t1.^2)
//disp(r2)
r3 = (-1)*intsplin(t1,theta1.*t1.^3)
//disp(r3)

//find Z
j=(r2^2)/4-((1/6)*r1*r3);
k=((r1*r2)/2)-((r0*r3)/6);
l=(r1^2)-((1/2)*r0*r2);
m=(k^2)-(4*j*l);
x1=(-k+sqrt(m))/(2*j);
x2=(-k-sqrt(m))/(2*j);
printf("Z = %f\n" ,x1)
printf("Z = %f\n" ,x2)

//find p
p1=-(r0+(x2)*r1)/(r1+(x2/2)*r2);
printf("p = %f\n" ,p1)

p2=-(r0+(x1)*r1)/(r1+(x1/2)*r2);
printf("p = %f\n" ,p2)

//find alpha
alpha=-(p1^2)*(r0+(x2*r1));
printf("alpha = %f\n" ,alpha)

t_d = -(1/(p1-x2))*log(((x2^2)*(alpha+(p1*(p1-x2)*r0)))/(alpha*(p1^2)))
printf("td = %f\n" ,t_d)


