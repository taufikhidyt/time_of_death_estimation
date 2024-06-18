clc;
clear;
 
//Ts=15;
Z=0.02;
p=0.8;
Ts = 35
td=-3.5;
Td=38;
t=0:2.5:10;
T=Ts+(Td-Ts)*((p/(p-Z))*exp(-Z*(t-td))-(Z/(p-Z)*exp(-p*(t-td))));
theta = T-Ts
disp(theta);

//find initial guess Z0 and p0
Z0=(1/(t($)-t($-1)))*log(theta($-1)/theta($));
disp(Z0);
a0=theta($)*exp(Z0*t($))-theta(1);
//disp(a0);
p0=(theta(2)-theta(1)+(Z0*t(2)*theta($)*exp(Z0*t($))))/(t(2)*(theta($)*exp(Z0*t($))-theta(1)));
disp(p0);

//Define the objective function
function E = objFunction(x)
    n = length(t);      // Get the length of the time array
    u = zeros(1, n);    // Initialize the vector u with zeros
    v = zeros(1, n);    // Initialize the vector v with zeros
    
    // Loop over the elements and calculate v
    for i = 1:n
        u(i) = exp(-x(2)*t(i)) - exp(-x(1)*t(i));
        v(i) = theta(1)*exp(-x(1)*t(i)) - theta(i)
    end

    // Return the vector
    u;
    v;
    
    // create the objective function
    n1 = norm(v);
    n2 = norm(u);
    dot = sum(u .* v);
    E = n1^2 - ((dot)^2)/((n2)^2) //objective function
endfunction

x0 = [p0; Z0]; // Initial guess for the solution
xopt = fminsearch(objFunction, x0); // Perform the minimization
disp(xopt)
p1 = xopt(1)
Z1 = xopt(2)
//disp(Z1)

// find new a
a1 = theta($)*exp(Z1*t($))-theta(1);
disp(a1)

// find new td
t_d = (1/(p1-Z1))*log((a1*p1)/((theta(1)+a1)*Z1))
disp(t_d)

