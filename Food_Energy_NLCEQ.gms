*----------------------------------------------------------------------
* This program solves the dynamic stochastic model of food and clean energy using NLCEQ
*
* Author: Yongyang Cai, University of Chicago and Hoover Institution
*
* If using material from this code, the user should cite the following paper:
* Cai, Y., K.L. Judd, and J. Steinbuks (2016). A nonlinear certainty equivalent* approximation method for dynamic stochastic problems. Quantitative Economics (forthcoming).
*----------------------------------------------------------------------

Set j land types: "1" = Food "2" = Biofuels /1*2/;

scalars
beta0 discount rate /0.95/
alpha biofuels cost share /0.5/
delta polution decay /0.001/
mu oil emissions factor relative to biofuels /0.25/
rhoe energy CES function parameter /0.5/
gamma intertemporal elasticity of substitution /0.5/
Be relative weight for energy production /0.5/
Bo relative weight for other goods /0.5/
BM relative weight for pollution disutility / 1 /
ae technology parameter for energy /1/
psio1 unit cost of oil extraction {psi^{o}_1} /0.4/
psio2 cost of oil extraction parameter {psi^{o}_2} /1/
eta parameter of disutility of pollution {eta} /4/
X0 Oil reserves in time 1 /1/
M0 Accumulated pollution in time 1 /1/
Lbar Total land amount /1/
target Pollution target /1.06/
gammahat utility parameter
epsilon a small number / 0.000001 /
;

gammahat = 1-(1/gamma);

Parameters
c(j) unit costs of food and biofuels {psi^{f} and psi^{Be}}
;

c("1") = 0.3;
c("2") = 0.5;

set n index of yields of food /1*2/;
alias(n, n2, n3, n4);

parameter foodYield(n);
foodYield('1') = 1;
foodYield('2') = 0.9;

parameter p21 probability of tipping in one year / 0.0034 /;

parameter tranProbs(n,n2);
tranProbs(n,n2) = 0;
tranProbs('1','1') = 1-p21;
tranProbs('2','1') = p21;
tranProbs('1','2') = 0;
tranProbs('2','2') = 1;

*----------------------------------------------------------------------
* Step 1: Define the transformed deterministic model

set t time /1*201/;

parameters
beta(t) discount rate
foodYieldPath(t)
probs(t,n)
;
beta(t) = beta0**(ord(t) - 1);

probs(t,n) = 0;
probs('1','1') = 1;
foodYieldPath(t) = 1;

loop(t$(ord(t)<card(t)),
probs(t+1,n) = sum(n2, tranProbs(n,n2) * probs(t,n2));
foodYieldPath(t+1) = sum(n2, probs(t+1,n2)*foodYield(n2));
);



Variables
obj objective criterion
u(t) utility
termU
;

Positive variables
M(t) pollution stock
X(t) oil stock
dx(t) oil extracted
ye(t) energy
yeTerm
yo(t) other goods
yoTerm
L(j,t) land allocation
;

Equations
tland(t) land constraint
energy(t) energy production function
energyTerm
othergoods(t) other goods
othergoodsTerm
oilextr(t) oil extraction function
pollution(t) pollution accumulation function
util(t) utility function
utilTerm
obj1 objective function
;

tland(t) .. sum(j, L(j,t)) =l= Lbar;

energy(t)$(ord(t) lt card(t)).. 
ye(t) =e= ae*( alpha * ( L("2",t) )**rhoe + (1-alpha) * (dx(t)+epsilon)**rhoe )**(1/rhoe);

energyTerm.. 
yeTerm =e= ae*( alpha * ( sum(t$(ord(t)=card(t)),L("2",t)) )**rhoe + (1-alpha) * (sum(t$(ord(t)=card(t)),0.01*X(t))+epsilon)**rhoe )**(1/rhoe);

oilextr(t)$(ord(t) lt card(t)).. 
X(t+1) =e= X(t) - dx(t);

pollution(t)$(ord(t) lt card(t)).. 
M(t+1) =e= (1-delta)*M(t) + mu*dx(t);

othergoods(t)$(ord(t) lt card(t)).. 
yo(t) =e= 1 - sum(j, c(j)*L(j,t)) - psio1*dx(t)/((X(t) + epsilon)**psio2);

othergoodsTerm.. 
yoTerm =e= 1 - sum(j, c(j)*sum(t$(ord(t)=card(t)),L(j,t))) - psio1*0.01*sum(t$(ord(t)=card(t)),X(t)/((X(t) + epsilon)**psio2));

util(t)$(ord(t) lt card(t)) ..
u(t) =e= ((foodYieldPath(t)*L("1",t))**gammahat)/gammahat
  + Be * (ye(t)**gammahat)/gammahat + Bo * (yo(t)**gammahat)/gammahat 
  - BM * M(t)**eta;

utilTerm ..
termU =e= ((sum(t$(ord(t)=card(t)),foodYieldPath(t)*L("1",t)))**gammahat)/gammahat
  + Be * (yeTerm**gammahat)/gammahat + Bo * (yoTerm**gammahat)/gammahat 
  - BM * (sum(t$(ord(t)=card(t)),M(t)))**eta;

obj1 .. obj =e= sum(t$(ord(t) lt card(t)), beta(t)*(u(t))) + sum(t$(ord(t)=card(t)), beta(t))*termU/(1-beta0);


* bound constraints

L.up(j,t) = Lbar;
M.up(t) = target;
L.lo(j,t) = epsilon;
yo.lo(t) = epsilon;
ye.lo(t) = epsilon;
M.lo(t) = epsilon;
yoTerm.lo = epsilon;
yeTerm.lo = epsilon;

* Initial Guess

L.l("1",t) = 0.5;
L.l("2",t) = Lbar - L.l("1",t);
X.l('1') = X0;
dx.l(t) = X.l(t)/card(t);
loop (t,
  X.l(t+1) = X.l(t) - dx.l(t);
);
ye.l(t) = ae*( alpha * ( L.l("2",t) )**rhoe + (1-alpha) * (dx.l(t)+epsilon)**rhoe )**(1/rhoe);
yeTerm.l = ae*( alpha * ( sum(t$(ord(t)=card(t)),L.l("2",t)) )**rhoe + (1-alpha) * (sum(t$(ord(t)=card(t)),0.01*X.l(t))+epsilon)**rhoe )**(1/rhoe);
yo.l(t) = 1 - sum(j, c(j)*L.l(j,t)) - psio1*dx.l(t)/((X.l(t) + epsilon)**psio2);
yoTerm.l = 1 - sum(j, c(j)*sum(t$(ord(t)=card(t)),L.l(j,t))) - psio1*0.01*sum(t$(ord(t)=card(t)),X.l(t)/((X.l(t) + epsilon)**psio2));
M.l('1') = M0;
loop (t,
  M.l(t+1) = (1-delta)*M.l(t) + mu*dx.l(t);
);
u.l(t)$(ord(t) lt card(t)) =((foodYieldPath(t)*L.l("1",t))**gammahat)/gammahat
  + Be * (ye.l(t)**gammahat)/gammahat + Bo * (yo.l(t)**gammahat)/gammahat
  - BM * M.l(t)**eta;
termU.l = ((sum(t$(ord(t)=card(t)),foodYieldPath(t)*L.l("1",t)))**gammahat)/gammahat
  + Be * (yeTerm.l**gammahat)/gammahat + Bo * (yoTerm.l**gammahat)/gammahat 
  - BM * (sum(t$(ord(t)=card(t)),M.l(t)))**eta;
obj.l = sum(t$(ord(t) lt card(t)), beta(t)*(u.l(t))) + sum(t$(ord(t)=card(t)), beta(t))*termU.l/(1-beta0);


model land / all /;

Option limrow=0,limcol=0,solprint=off;
option decimals = 7;
OPTION ITERLIM = 500000;
OPTION RESLIM = 500000;
OPTION DOMLIM = 1000;
OPTION NLP= conopt;


*----------------------------------------------------------------------
* Step 2: Optimization step

set v1 grid space / 1*21/;
alias(v1,v2,d1,d2);

Parameters
Xmin    lower bound of X  / 0.01 /
Xmax    upper bound of X  / 1 /
pmin    lower bound of M  / 1 /
pmax    upper bound of M  / 1.06 /

Xinit(v1) Grid: Oil stock
Minit(v2) Grid: Pollution stock
;

* Chebyshev parameters

parameters zx(v1)
        zp(v2)
        extx
        extp
        extxmin
        extpmin
        extxmax
        extpmax
        dzx
        dzp
;

zx(v1) = - cos((2*ord(v1)-1)/(2*card(v1))*pi);
zp(v2) = - cos((2*ord(v2)-1)/(2*card(v2))*pi);
extx = -(zx('1')+1)*(Xmax-Xmin)/2/zx('1');
extp = -(zp('1')+1)*(pmax-pmin)/2/zp('1');
extxmin = Xmin - extx;
extpmin = pmin - extp;
extxmax = Xmax + extx;
extpmax = pmax + extp;
dzx = 2/(extxmax-extxmin);
dzp = 2/(extpmax-extpmin);


Xinit(v1) = extxmin + (1+zx(v1))/dzx;
Minit(v2) = extpmin + (1+zp(v2))/dzp;


Parameters
valf(v1,v2,n) value function
Lfsol(v1,v2,n)
Lbsol(v1,v2,n)
dxsol(v1,v2,n)
;

loop ((v1,v2,n),

  X.fx(t)$(ord(t) eq 1) = Xinit(v1);
  M.fx(t)$(ord(t) eq 1) = Minit(v2);
  foodYieldPath('1') = foodYield(n);

  probs(t,n3) = 0;
  probs('1',n) = 1;
  loop(t$(ord(t)<card(t)),
    probs(t+1,n3) = sum(n2, tranProbs(n3,n2) * probs(t,n2));
    foodYieldPath(t+1) = sum(n3, probs(t+1,n3)*foodYield(n3));
  );

  solve land using nlp maximizing obj;

  valf(v1,v2,n) = obj.l;
  Lfsol(v1,v2,n) = L.l('1','1');
  Lbsol(v1,v2,n) = L.l('2','1');
  dxsol(v1,v2,n) = dx.l('1');
);


*----------------------------------------------------------------------
* Step 3: Approximation step


parameters TTx(v1,d1)
        TTp(v2,d2)
;
TTx(v1,'1') = 1;
TTx(v1,'2') = zx(v1);
loop(d1$(ord(d1)>=2 and ord(d1)<card(d1)),
  TTx(v1,d1+1) = 2*zx(v1)*TTx(v1,d1) - TTx(v1,d1-1);
);
TTp(v2,'1') = 1;
TTp(v2,'2') = zp(v2);
loop(d2$(ord(d2)>=2 and ord(d2)<card(d2)),
  TTp(v2,d2+1) = 2*zp(v2)*TTp(v2,d2) - TTp(v2,d2-1);
);

parameter coefs(d1,d2,n) approximation coefficients;
parameter coefsLf(d1,d2,n) approximation coefficients;
parameter coefsLb(d1,d2,n) approximation coefficients;
parameter coefsdx(d1,d2,n) approximation coefficients;

coefs(d1,d2,n) = 0;
coefs(d1,d2,n)$(ord(d1)+ord(d2)<=card(d1)+1) =
    4*sum((v1,v2), TTx(v1,d1)*TTp(v2,d2)*valf(v1,v2,n)) / (card(v1)*card(v2));
coefs('1',d2,n) = coefs('1',d2,n)/2;
coefs(d1,'1',n) = coefs(d1,'1',n)/2;

coefsLf(d1,d2,n) = 0;
coefsLf(d1,d2,n)$(ord(d1)+ord(d2)<=card(d1)+1) =
    4*sum((v1,v2), TTx(v1,d1)*TTp(v2,d2)*Lfsol(v1,v2,n)) / (card(v1)*card(v2));
coefsLf('1',d2,n)= coefsLf('1',d2,n)/2;
coefsLf(d1,'1',n) = coefsLf(d1,'1',n)/2;

coefsLb(d1,d2,n) = 0;
coefsLb(d1,d2,n)$(ord(d1)+ord(d2)<=card(d1)+1) =
    4*sum((v1,v2), TTx(v1,d1)*TTp(v2,d2)*Lbsol(v1,v2,n)) / (card(v1)*card(v2));
coefsLb('1',d2,n)= coefsLb('1',d2,n)/2;
coefsLb(d1,'1',n) = coefsLb(d1,'1',n)/2;

coefsdx(d1,d2,n) = 0;
coefsdx(d1,d2,n)$(ord(d1)+ord(d2)<=card(d1)+1) =
    4*sum((v1,v2), TTx(v1,d1)*TTp(v2,d2)*dxsol(v1,v2,n)) / (card(v1)*card(v2));
coefsdx('1',d2,n)= coefsdx('1',d2,n)/2;
coefsdx(d1,'1',n) = coefsdx(d1,'1',n)/2;

display coefs;

*----------------------------------------------------------------------
* Step 4: Output solutions for error checking

set k1 /1*100/;
alias(k1, k2);

Parameters
valf1(k1,k2,n) value function
Lfsol1(k1,k2,n)
Lbsol1(k1,k2,n)
dxsol1(k1,k2,n)
X1(k1) Grid: Oil stock
p1(k2) Grid: Pollution stock
;

X1(k1) = Xmin + (Xmax-Xmin)*(ord(k1)-1) / (card(k1)-1);
p1(k2) = pmin + (pmax-pmin)*(ord(k2)-1) / (card(k2)-1);
valf1(k1,k2,n) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefs(d1,d2,n)*cos((ord(d1)-1)*arccos(dzx*(X1(k1)-extxmin)-1)) * cos((ord(d2)-1)*arccos(dzp*(p1(k2)-extpmin)-1)) );
Lfsol1(k1,k2,n) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLf(d1,d2,n)*cos((ord(d1)-1)*arccos(dzx*(X1(k1)-extxmin)-1)) * cos((ord(d2)-1)*arccos(dzp*(p1(k2)-extpmin)-1)) );
Lbsol1(k1,k2,n) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLb(d1,d2,n)*cos((ord(d1)-1)*arccos(dzx*(X1(k1)-extxmin)-1)) * cos((ord(d2)-1)*arccos(dzp*(p1(k2)-extpmin)-1)) );
dxsol1(k1,k2,n) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsdx(d1,d2,n)*cos((ord(d1)-1)*arccos(dzx*(X1(k1)-extxmin)-1)) * cos((ord(d2)-1)*arccos(dzp*(p1(k2)-extpmin)-1)) );


File modeloutput /sol_cf.csv/;
modeloutput.nw=18;
modeloutput.nr=2;
modeloutput.nz=1e-15;

Put modeloutput;
Put "FoodYield, OilStock, Pollution, LandFood, LandBiofuel, OilExtraction, V" /;

modeloutput.pc=5;
*modeloutput.pw=4000;

loop((n,k1,k2),
  put foodYield(n):14:6;
  put X1(k1):14:6;
  put p1(k2):14:6;
  put Lfsol1(k1,k2,n):14:6;
  put Lbsol1(k1,k2,n):14:6;
  put dxsol1(k1,k2,n):14:6;
  put valf1(k1,k2,n):14:6;
  put /;
);


*----------------------------------------------------------------------
* Step 5: Output solutions over a path for plotting figures

set k4 /1*201/;

Parameters
aaM(k4,n)
aaX(k4,n)
aavalf(k4,n) value function
aaLfsol(k4,n)
aaLbsol(k4,n)
aadxsol(k4,n)
;

aaM('1',n) = M0;
aaX('1',n) = X0;
loop((k4,n)$(ord(k4)<card(k4)),
  aavalf(k4,n) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefs(d1,d2,n)*cos((ord(d1)-1)*arccos(dzx*(aaX(k4,n)-extxmin)-1)) * cos((ord(d2)-1)*arccos(dzp*(aaM(k4,n)-extpmin)-1)) );
  aaLfsol(k4,n) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLf(d1,d2,n)*cos((ord(d1)-1)*arccos(dzx*(aaX(k4,n)-extxmin)-1)) * cos((ord(d2)-1)*arccos(dzp*(aaM(k4,n)-extpmin)-1)) );
  aaLbsol(k4,n) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsLb(d1,d2,n)*cos((ord(d1)-1)*arccos(dzx*(aaX(k4,n)-extxmin)-1)) * cos((ord(d2)-1)*arccos(dzp*(aaM(k4,n)-extpmin)-1)) );
  aadxsol(k4,n) = sum( (d1,d2)$(ord(d1)+ord(d2)<=card(d1)+1),
        coefsdx(d1,d2,n)*cos((ord(d1)-1)*arccos(dzx*(aaX(k4,n)-extxmin)-1)) * cos((ord(d2)-1)*arccos(dzp*(aaM(k4,n)-extpmin)-1)) );

* it has some small approximation errors so that aaM could be slightly bigger than its upper bound,
* and then in the next period the arccos function in the approximation will not be well defined,
* so here we adjust it to not exceed its upper bound in order to proceed in next periods
  aaM(k4+1,n) = min(target, (1-delta)*aaM(k4,n) + mu*aadxsol(k4,n));

  aaX(k4+1,n) = aaX(k4,n) - aadxsol(k4,n);
);

File modeloutput3 /sol_cf_path.csv/;
modeloutput3.nw=18;
modeloutput3.nr=2;
modeloutput3.nz=1e-15;

Put modeloutput3;
Put "FoodYield, OilStock, Pollution, LandFood, LandBiofuel, OilExtraction, V" /;

modeloutput3.pc=5;
*modeloutput3.pw=4000;

loop((n,k4),
  put foodYield(n):14:6;
  put aaX(k4,n):14:6;
  put aaM(k4,n):14:6;
  put aaLfsol(k4,n):14:6;
  put aaLbsol(k4,n):14:6;
  put aadxsol(k4,n):14:6;
  put aavalf(k4,n):14:6;
  put /;
);

