*----------------------------------------------------------------------
* This program solves the RBC model with a constraint on investment
*
* Author: Yongyang Cai, University of Chicago and Hoover Institution
*
* If using material from this code, the user should cite the following paper:
* Cai, Y., K.L. Judd, and J. Steinbuks (2015). A nonlinear certainty equivalent* approximation method for dynamic stochastic problems. Quantitative Economics (forthcoming).
*----------------------------------------------------------------------

scalars
beta 	discount rate / 0.96 /
alpha 	capital share / 0.33 /
delta 	capital stock depreciation / 0.1 /
rho 	persistence parameter for TFP growth / 0.9 /
sigma   standard deviation of technology innovation / 0.013 /
phi 	threshold of investment constraint / 0.975 /
gamma 	relative risk aversion / 2 /
kss	steady state
Iss	investment in the steady state
;

kss = ((1/beta-1+delta)/alpha)**(1/(alpha-1));
Iss = delta*kss;

display Iss;

*************************************************
* Step 1: Define the transformed deterministic problem

Sets t time /1*101/;

parameters
a(t)		productivity path
betas(t)	discounts
;
betas(t) = beta**(ord(t)-1);
a(t) = 1;

Variables
obj objective criterion
u(t) utility
i(t) investment
;

Positive variables
k(t) capital stock
c(t) consumption
;

Equations

obj1 Objective function
util(t) Utility function
cap(t) Law of Motion for Capital Stock
bc(t) budget constraint
investBound(t) investment constraint
;

cap(t)$(ord(t) lt card(t)) .. k(t+1) =e= (1-delta)*k(t) + i(t);
bc(t)$(ord(t) lt card(t)) .. c(t) + i(t) =e= a(t)*k(t)**alpha;
util(t)$(ord(t) lt card(t)) .. u(t) =e= (c(t)**(1-gamma)-1)/(1-gamma);
investBound(t)$(ord(t) lt card(t)) .. i(t) =g= phi*Iss;
obj1 .. obj =e= sum(t$(ord(t) lt card(t)), betas(t)*u(t)) + 
  sum(t$(ord(t)=card(t)), betas(t)/(1-beta) * ( (0.7*a(t)*k(t)**alpha)**(1-gamma)-1 )/(1-gamma) );

* Constraints

k.lo(t) = 0.01;
k.up(t) = 100;
c.lo(t) = 0.001;
c.up(t) = 100;

* Initial Guess

i.l(t) = Iss;
k.l(t) = kss;
c.l(t) = a(t)*k.l(t)**alpha - i.l(t);
u.l(t) =(c.l(t)**(1-gamma)-1)/(1-gamma);
obj.l = sum(t$(ord(t) lt card(t)), betas(t)*u.l(t)) + 
  sum(t$(ord(t)=card(t)), betas(t)/(1-beta) * ( (0.7*a(t)*k.l(t)**alpha)**(1-gamma)-1 )/(1-gamma) );


Option limrow=0,limcol=0,solprint=off;
OPTION ITERLIM = 500000;
OPTION RESLIM = 500000;
OPTION DOMLIM = 1000;
OPTION NLP= conopt;
option decimals = 7;

model RBC /all/;


*******************************************************
* Step 2: Optimization step

set v1 grid space /1*21/;
alias(v1,v2,d1,d2);

* approximation domain
parameters
amin    lower bound of A  / 0.5 /
amax    upper bound of A  / 1.5 /
kmin    lower bound of k  
kmax    upper bound of k  
;

kmin = 0.5*kss;
kmax = 1.5*kss;

* Chebyshev parameters

parameters zk(v1)
        za(v2)
        extk
        exta
        extkmin
        extamin
        extkmax
        extamax
        dzk
        dza
;

zk(v1) = - cos((2*ord(v1)-1)/(2*card(v1))*pi);
za(v2) = - cos((2*ord(v2)-1)/(2*card(v2))*pi);
extk = -(zk('1')+1)*(kmax-kmin)/2/zk('1');
exta = -(za('1')+1)*(amax-amin)/2/za('1');
extkmin = kmin - extk;
extamin = amin - exta;
extkmax = kmax + extk;
extamax = amax + exta;
dzk = 2/(extkmax-extkmin);
dza = 2/(extamax-extamin);

Parameters
kinit(v1) Grid: Capital Stock 
ainit(v2) Grid: Productivity  
;

kinit(v1) = extkmin + (1+zk(v1))/dzk;
ainit(v2) = extamin + (1+za(v2))/dza;

Parameters
csol(v1,v2) solution for c
lambda(v1,v2) multiplier for the investment constraint
;

loop ((v1,v2),
  k.fx('1') = kinit(v1);
  a('1') = ainit(v2);
  loop (t$(ord(t) lt card(t)),
    a(t+1) = a(t)**rho;
  );

  solve RBC using nlp maximizing obj;

  lambda(v1,v2) = abs(investBound.m('1'));
  csol(v1,v2) = c.l('1');
);


************************************************
* Step 3: Approximation step

parameters TTk(v1,d1)
        TTa(v2,d2)
;
TTk(v1,'1') = 1;
TTk(v1,'2') = zk(v1);
loop(d1$(ord(d1)>=2 and ord(d1)<card(d1)),
  TTk(v1,d1+1) = 2*zk(v1)*TTk(v1,d1) - TTk(v1,d1-1);
);
TTa(v2,'1') = 1;
TTa(v2,'2') = za(v2);
loop(d2$(ord(d2)>=2 and ord(d2)<card(d2)),
  TTa(v2,d2+1) = 2*za(v2)*TTa(v2,d2) - TTa(v2,d2-1);
);

parameter coefsLam(d1,d2) approximation coefficients;
parameter coefsC(d1,d2) approximation coefficients;

coefsLam(d1,d2) = 0;
coefsLam(d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1) =
    4*sum((v1,v2), TTk(v1,d1)*TTa(v2,d2)*lambda(v1,v2)) / (card(v1)*card(v2));
coefsLam('1',d2) = coefsLam('1',d2)/2;
coefsLam(d1,'1') = coefsLam(d1,'1')/2;

coefsC(d1,d2) = 0;
coefsC(d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1) =
    4*sum((v1,v2), TTk(v1,d1)*TTa(v2,d2)*csol(v1,v2)) / (card(v1)*card(v2));
coefsC('1',d2)= coefsC('1',d2)/2;
coefsC(d1,'1') = coefsC(d1,'1')/2;


display coefsC, coefsLam;


*******************************************************
* Step 4: Error Checking

* Compute Approximation Errors

set v3 / 1*1000 /;
Parameters
akinit(v3) Capital Stock 
aainit(v3) Productivity  
;

akinit(v3) = uniform(kmin,kmax);
aainit(v3) = uniform(amin,amax);

Parameters
acsol(v3) solution for c
alambda(v3) multiplier for the investment constraint
;

loop (v3,
  k.fx('1') = akinit(v3);
  a('1') = aainit(v3);

  loop (t$(ord(t) lt card(t)),
    a(t+1) = a(t)**rho;
  );

  solve RBC using nlp maximizing obj;

  alambda(v3) = abs(investBound.m('1'));
  acsol(v3) = c.l('1');
);

parameters
errsLam(v3)
errsC(v3)
errLamLinf
errLamL1
errCLinf
errCL1
;

errsLam(v3) = abs(alambda(v3) - sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLam(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(akinit(v3)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aainit(v3)-extamin)-1)) ) ) / (1+abs(alambda(v3)));
errsC(v3) = abs(acsol(v3) - sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsC(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(akinit(v3)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aainit(v3)-extamin)-1)) ) ) / (1+abs(acsol(v3)));
errLamLinf = smax(v3, errsLam(v3));
errLamL1 = sum(v3, errsLam(v3)) / card(v3);
errCLinf = smax(v3, errsC(v3));
errCL1 = sum(v3, errsC(v3)) / card(v3);

display errLamLinf, errLamL1, errCLinf, errCL1;


*******************************************************
* Compute global errors

set k1 /1*10000/;
set k3 index of Gauss-Hermite quadrature rule /1*15/; 

Parameters
lambda1(k1) 
lambda1p(k1,k3) 
kpsol1(k1)
csol1(k1)
csol1p(k1,k3)
Isol1(k1)
kmin1
kmax1
amin1
amax1
kk(k1) capital
aa(k1) productivity 
aap(k1,k3)
integrand(k1,k3)
;


parameter GHNodes(k3) Gaussian Hermite quadrature nodes;
GHNodes('1') =  -4.499990707309392*sqrt(2.0);
GHNodes('2') =  -3.669950373404453*sqrt(2.0);
GHNodes('3') =  -2.967166927905603*sqrt(2.0);
GHNodes('4') =  -2.325732486173858*sqrt(2.0);
GHNodes('5') =  -1.719992575186489*sqrt(2.0);
GHNodes('6') =  -1.136115585210921*sqrt(2.0);
GHNodes('7') =  -0.5650695832555757*sqrt(2.0);
GHNodes('8') =  0;
GHNodes('9') =  -GHNodes('7') ;
GHNodes('10') =  -GHNodes('6') ;
GHNodes('11') =  -GHNodes('5') ;
GHNodes('12') =  -GHNodes('4') ;
GHNodes('13') =  -GHNodes('3') ;
GHNodes('14') =  -GHNodes('2') ;
GHNodes('15') =  -GHNodes('1') ;

parameter GHWeights(k3) Gaussian Hermite quadrature weights;
GHWeights('1') =  1.522475804253517e-9/sqrt(pi);
GHWeights('2') =  1.059115547711067e-6/sqrt(pi);
GHWeights('3') =  0.0001000044412324999/sqrt(pi);
GHWeights('4') =  0.002778068842912776/sqrt(pi);
GHWeights('5') =  0.03078003387254608/sqrt(pi);
GHWeights('6') =  0.1584889157959357/sqrt(pi);
GHWeights('7') =  0.4120286874988986/sqrt(pi);
GHWeights('8') =  0.5641003087264175/sqrt(pi);
GHWeights('9') =  GHWeights('7');
GHWeights('10') =  GHWeights('6');
GHWeights('11') =  GHWeights('5');
GHWeights('12') =  GHWeights('4');
GHWeights('13') =  GHWeights('3');
GHWeights('14') =  GHWeights('2');
GHWeights('15') =  GHWeights('1');

kmin1 = 0.7*kss;
kmax1 = 1.3*kss;
amin1 = 0.7;
amax1 = 1.3;

kk(k1) = Uniform(kmin1, kmax1);
aa(k1) = Uniform(amin1, amax1);
aap(k1,k3) = aa(k1)**rho * exp(sigma*GHNodes(k3));
lambda1(k1) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLam(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kk(k1)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aa(k1)-extamin)-1)) );
csol1(k1) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsC(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kk(k1)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aa(k1)-extamin)-1)) );
Isol1(k1) = aa(k1)*kk(k1)**alpha - csol1(k1);
kpsol1(k1) = (1-delta)*kk(k1) + Isol1(k1);
csol1p(k1,k3) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsC(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kpsol1(k1)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aap(k1,k3)-extamin)-1)) );
lambda1p(k1,k3) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLam(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kpsol1(k1)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aap(k1,k3)-extamin)-1)) );

integrand(k1,k3) = csol1p(k1,k3)**(-gamma) * ( 1-delta+alpha*aap(k1,k3)*kpsol1(k1)**(alpha-1) ) - (1-delta)*lambda1p(k1,k3);

parameters
errs1(k1)
errs2(k1)
errs3(k1)
err1Linf
err2Linf
err3Linf
err1L1
err2L1
err3L1
;

errs1(k1) = abs((lambda1(k1) + beta*sum(k3, GHWeights(k3)*integrand(k1,k3))) / csol1(k1)**(-gamma) - 1);
err1Linf = smax(k1, errs1(k1));
err1L1 = sum(k1, errs1(k1)) / card(k1);

errs2(k1) = abs(lambda1(k1) * (Isol1(k1)/(phi*Iss)-1));
err2Linf = smax(k1, errs2(k1));
err2L1 = sum(k1, errs2(k1)) / card(k1);

errs3(k1) = max(0, 1-Isol1(k1)/(phi*Iss));
err3Linf = smax(k1, errs3(k1));
err3L1 = sum(k1, errs3(k1)) / card(k1);

display err1Linf, err1L1, err2Linf, err2L1, err3Linf, err3L1;

*****************************************************
* Compute weighted errors with ergodic points

parameters
werrs1(k1)
werrs2(k1)
werrs3(k1)
werr1Linf
werr2Linf
werr3Linf
werr1L1
werr2L1
werr3L1
;

kk('1') = kss;
aa('1') = 1;
loop(k1$(ord(k1)<card(k1)),
  csol1(k1) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsC(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kk(k1)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aa(k1)-extamin)-1)) );
  Isol1(k1) = aa(k1)*kk(k1)**alpha - csol1(k1);
  aa(k1+1) = aa(k1)**rho * exp(sigma*normal(0,1));
  kk(k1+1) = kk(k1)*(1-delta) + Isol1(k1);
);

aap(k1,k3) = aa(k1)**rho * exp(sigma*GHNodes(k3));
lambda1(k1) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLam(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kk(k1)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aa(k1)-extamin)-1)) );
csol1(k1) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsC(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kk(k1)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aa(k1)-extamin)-1)) );
Isol1(k1) = aa(k1)*kk(k1)**alpha - csol1(k1);
kpsol1(k1) = (1-delta)*kk(k1) + Isol1(k1);
csol1p(k1,k3) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsC(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kpsol1(k1)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aap(k1,k3)-extamin)-1)) );
lambda1p(k1,k3) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLam(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kpsol1(k1)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aap(k1,k3)-extamin)-1)) );

integrand(k1,k3) = csol1p(k1,k3)**(-gamma) * ( 1-delta+alpha*aap(k1,k3)*kpsol1(k1)**(alpha-1) ) - (1-delta)*lambda1p(k1,k3);


werrs1(k1) = abs((lambda1(k1) + beta*sum(k3, GHWeights(k3)*integrand(k1,k3))) / csol1(k1)**(-gamma) - 1);
werr1Linf = smax(k1, werrs1(k1));
werr1L1 = sum(k1, werrs1(k1)) / card(k1);

werrs2(k1) = abs(lambda1(k1) * (Isol1(k1)/(phi*Iss)-1));
werr2Linf = smax(k1, werrs2(k1));
werr2L1 = sum(k1, werrs2(k1)) / card(k1);

werrs3(k1) = max(0, 1-Isol1(k1)/(phi*Iss));
werr3Linf = smax(k1, werrs3(k1));
werr3L1 = sum(k1, werrs3(k1)) / card(k1);

display werr1Linf, werr1L1, werr2Linf, werr2L1, werr3Linf, werr3L1;



*****************************************************
* Compute errors for plotting figures

set k4 /1*101/;
set k5 /1*3/; 

Parameters
kkgrid(k4)
csol3(k4,k5)
aa3(k5) productivity 
csol3(k4,k5)
aap3(k5,k3)
lambda3(k4,k5)
Isol3(k4,k5)
kpsol3(k4,k5)
csol3p(k4,k5,k3)
lambda3p(k4,k5,k3)
integrand3(k4,k5,k3)
aerrs1(k4,k5)
aerrs2(k4,k5)
aerrs3(k4,k5)
aerrs(k4,k5)
;

kkgrid(k4) = kmin1 + (ord(k4)-1)*(kmax1-kmin1)/(card(k4)-1);
aa3('1') = amin1;
aa3('2') = 1;
aa3('3') = amax1;

aap3(k5,k3) = aa3(k5)**rho * exp(sigma*GHNodes(k3));
lambda3(k4,k5) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLam(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kkgrid(k4)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aa3(k5)-extamin)-1)) );
csol3(k4,k5) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsC(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kkgrid(k4)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aa3(k5)-extamin)-1)) );

Isol3(k4,k5) = aa3(k5)*kkgrid(k4)**alpha - csol3(k4,k5);
kpsol3(k4,k5) = (1-delta)*kkgrid(k4) + Isol3(k4,k5);
csol3p(k4,k5,k3) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsC(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kpsol3(k4,k5)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aap3(k5,k3)-extamin)-1)) );
lambda3p(k4,k5,k3) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLam(d1,d2)*cos((ord(d1)-1)*arccos(dzk*(kpsol3(k4,k5)-extkmin)-1)) * cos((ord(d2)-1)*arccos(dza*(aap3(k5,k3)-extamin)-1)) );

integrand3(k4,k5,k3) = csol3p(k4,k5,k3)**(-gamma) * ( 1-delta+alpha*aap3(k5,k3)*kpsol3(k4,k5)**(alpha-1) ) - (1-delta)*lambda3p(k4,k5,k3);

aerrs1(k4,k5) = abs((lambda3(k4,k5) + beta*sum(k3, GHWeights(k3)*integrand3(k4,k5,k3))) / csol3(k4,k5)**(-gamma) - 1);
aerrs2(k4,k5) = abs(lambda3(k4,k5) * (Isol3(k4,k5)/(phi*Iss)-1));
aerrs3(k4,k5) = max(0, 1-Isol3(k4,k5)/(phi*Iss));

aerrs(k4,k5) = max(aerrs1(k4,k5), max(aerrs2(k4,k5), aerrs3(k4,k5)));

File outputerr /growth_bind_cf_errs.csv/;
outputerr.nw=18;
outputerr.nr=2;
outputerr.nz=1e-15;
Put outputerr;
outputerr.pc=5;

loop(k4,
  put kkgrid(k4):14:6;
  loop(k5,
    put aerrs(k4,k5):14:6;
  );
  put /;
);
