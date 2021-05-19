
*Define scalars
scalars mu      gravitanional constant*mass      /1.327e20/
        r0      initial radius of orbit          /1.496e11/
        m0      initial mass of satellite        /4.53e3/
        mdot    massflow                         /6.76e-5/
        mfuel   mass of fuel
        tf      time in final orbit              /1.668e7/
        Thrust                                   /3.77/
        c3;


mfuel = 0.2*m0;
c3 =(tf*sqrt(mu))/sqrt(r0**3);


*Define time

Set k /0*500/;

Scalar deltat;
deltat = 1/(card(k));


*Define variables

variables   z
            x1(k)
            x2(k)
            x3(k)
            x1dot(k) ratio between initial and final tangential velocities
            x2dot(k)
            x3dot(k)
            u(k)    angle
            m(k)    mass
            p(k)    proportionality constant
            c1(k)
            c2(k);
            
positive variable p(k);
            
Parameter t(k) Are we using this?;
t(k) = deltat * (ord(k)-1); 


*Set boundaries
x1.lo(k) = 1;
m.lo(k) = 0.1;
u.lo(k) = -pi/2;
u.up(k)= pi/2;


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
            mass(k)
            mfuelmax(k)
            pmax(k)
            minitial;
            
mfuelmax(k)..   m0-m(k) =L= mfuel;
mass(k+1)..     m(k+1)=E= m(k) - p(k)*deltat*mdot*tf;
pmax(k)..       p(k) =L= 1;

obj.. z =E= x1("500");

endcond1.. x1dot("500") =E= 0;
endcond2..  x3("500") =E= 1 / sqrt(x1("500"));

eq1(k+1).. x1(k+1) =E= x1(k) + deltat * x1dot(k+1);
eq2(k+1).. x2(k+1) =E= x2(k) + deltat * x2dot(k+1);
eq3(k+1).. x3(k+1) =E= x3(k) + deltat * x3dot(k+1);


diffeq1(k)..  x1dot(k)=E=c3**2 *x2(k);
diffeq2(k)..  x2dot(k)=E=((x3(k)**2)/x1(k)) - (1/(x1(k)**2)) + r0**2*Thrust*p(k)*sin(u(k)) / (mu*m(k));
diffeq3(k)..  x3dot(k)=E=-(mu*tf**2*x2(k)*x3(k))/(r0**3*x1(k))+(tf*sqrt(r0)*Thrust*p(k)*cos(u(k)))/(sqrt(mu)*m(k));

initial1.. x1("0") =E= 1;
initial2.. x2("0") =E= 0;
*if an object is orbiting a planet at a constant velocity then v = sqrt(GM/r), or v = sqrt(mu/r)
initial3.. x3("0") =E= 1;
minitial.. m("0") =E= m0;



model maxradius /ALL/;

SOLVE maxradius USING NLP maximizing z;

Display z.L, u.L, p.L;
