
*Define scalars
scalars mu 
        r0 
        m0 
        mdot 
        tf 
        Thrust
        c1
        c2
        c3;

mu = 1.327*10**20;
r0 = 1.496*10**11;
m0 = 4.53*10**3;
mdot = 6.76*10**(-5);
tf = 1.668*10**7;
Thrust = 3.77;
c1 = (Thrust*r0**2)/(mu*m0);
c2 =(mdot*sqrt(mu))/(Thrust*sqrt(r0));
c3 =(tf*sqrt(mu))/sqrt(r0**3);


*Define time

Set k /0*100/;

Scalar deltat;
deltat = 1/(card(k));

Parameter t(k);
t(k) = deltat * (ord(k)-1) ;

*Define variables

variables   z
            x1(k)
            x2(k)
            x3(k)
            x1dot(k)
            x2dot(k)
            x3dot(k)
            u(k)    angle
            m(k)    mass;

*Set boundaries
x1.lo(k) = 0.1;
*u.l(k) = pi;
x1dot.l(k) = 10000;
x2dot.l(k) = 10000;
x3dot.l(k) = 10000;

equations   obj,
            eq1(k)
            eq2(k)
            eq3(k)
            diffeq1(k)
            diffeq2(k)
            diffeq3(k)
            initial1
            initial2
            initial3
            endcond1
            endcond2
            mass(k);

obj.. z =E= x1("100");

endcond1.. x1dot("100") =E= 0;
endcond2..  x3("100") =E= 1 / sqrt(x1("100"));

eq1(k+1).. x1(k+1) =E= x1(k) + deltat * x1dot(k);
eq2(k+1).. x2(k+1) =E= x2(k) + deltat * x2dot(k);
eq3(k+1).. x3(k+1) =E= x3(k) + deltat * x3dot(k);


diffeq1(k)..  x1dot(k)=E=c3**2 *x2(k);
diffeq2(k)..  x2dot(k)=E=((x3(k)**2)/x1(k)) - (1/(x1(k)**2)) + sin(u(k)) / (1/c1 - c2*c3*t(k));
diffeq3(k)..  x3dot(k)=E=(-c3**2)*((x2(k)*x3(k))/x1(k)) + (c3*cos(u(k)) / (1/c1 - c2*c3*t(k)));

initial1.. x1("0") =E= 1;
initial2.. x2("0") =E= 0;
*if an object is orbiting a planet at a constant velocity then v = sqrt(GM/r), or v = sqrt(mu/r)
initial3.. x3("0") =E= 1;

mass(k).. m(k) =E= m0 - mdot * tf * t(k)


model maxradius /ALL/;

SOLVE maxradius USING NLP maximizing z;

Display z.L, u.L, x1.L;
