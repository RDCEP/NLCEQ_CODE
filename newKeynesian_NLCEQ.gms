*----------------------------------------------------------------------
* This program solves the New Keynesian DSGE model with zero lower bound using NLCEQ
*
* Author: Yongyang Cai, University of Chicago and Hoover Institution

* If using material from this code, the user should cite the following paper:
* Cai, Y., K.L. Judd, and J. Steinbuks (2016). A nonlinear certainty equivalent* approximation method for dynamic stochastic problems. Quantitative Economics (forthcoming).
*----------------------------------------------------------------------

scalars
eta	labor supply elasticity / 1 /
beta 	discount factor / 0.994 /
rho 	persistence of discount rate shock / 0.8 /
sigma   standard deviation of discount rate shock / 0.005 /
alpha 	elasticity of substitution across intermediate goods / 6 /
theta   Calvo parameter / 0.9 /
sg	ratio of govenment spending to output / 0.2 /
phi_pi  response to inflation (monetary policy rule) / 2.5 /
phi_y   response to output (monetary policy rule) / 0.25 /
pi_ss   steady state inflation / 1.005 /
;

parameters
rss	steady state interest rate
chi2ss
chi1ss
qss
zss
lss
css
yss	steady state output
vss	steady state of v
;

rss = pi_ss/beta-1;
chi2ss = 1/((1-sg)*(1-theta*beta*pi_ss**(alpha-1)));
qss = ((1-theta*pi_ss**(alpha-1))/(1-theta))**(1/(1-alpha));
chi1ss = chi2ss*qss*(alpha-1)/alpha;
vss = (1-theta)*qss**(-alpha) / (1-theta*pi_ss**alpha);
yss = (chi1ss*(1-theta*beta*pi_ss**alpha)/(vss**eta))**(1/(1+eta));
zss = rss;
lss = vss*yss;
css = (1-sg)*yss;

display rss, vss, yss, lss, css, chi1ss, chi2ss, qss;

*******************************************************
* Step 1: Define the transformed deterministic model 

set t time index /1*201/;

parameter betas(t) discount;
betas(t) = beta;

set iii /1*8/;

Variables
z(t)
obj
;

Positive variables
v(t)
y(t)
pi_t(t)
q(t)
chi1(t)
chi2(t)
epss(iii)
;

Equations
Objective
StateFun(t)
NotionalInterest(t)
QFun(t)
Chi1Fun(t)
Chi2Fun(t)
Chi1Chi2(t)
Euler(t)
TermChi1
TermChi2
TermY
Termv
;

Objective..
obj =e= sum(iii, epss(iii));

StateFun(t)$(ord(t)<card(t))..
v(t+1) =e= (1-theta)*power(q(t),-alpha) + theta*power(pi_t(t),alpha)*v(t);

NotionalInterest(t)..
z(t) =e= (1+rss)*(pi_t(t)/pi_ss)**phi_pi*(y(t)/yss)**phi_y - 1;

QFun(t)..
1 =e= power(q(t),alpha-1) * (1-theta*power(pi_t(t),alpha-1))/(1-theta);

Chi1Fun(t)$(ord(t)<card(t))..
chi1(t) =e= power(y(t),1+eta)*power(v(t+1),eta) + theta*betas(t+1)*power(pi_t(t+1),alpha)*chi1(t+1);

Chi2Fun(t)$(ord(t)<card(t))..
chi2(t) =e= 1/(1-sg) + theta*betas(t+1)*power(pi_t(t+1),alpha-1)*chi2(t+1);

Chi1Chi2(t)..
chi1(t) =e= q(t)*chi2(t)*(alpha-1)/alpha;

Euler(t)$(ord(t)<card(t))..
pi_t(t+1)*y(t+1) =e= betas(t+1)*(1+max(0,z(t)))*y(t);

TermChi1..
chi1('201')/chi1ss-1 =e= epss('1') - epss('2');

TermChi2..
chi2('201')/chi2ss-1 =e= epss('3') - epss('4');

TermY..
y('201')/yss-1 =e= epss('5') - epss('6');

Termv..
v('201')/vss-1 =e= epss('7') - epss('8');

* bounds
y.lo(t) = 0.001;
v.lo(t) = 0.001;
pi_t.lo(t) = 0.001;
q.lo(t) = 0.001;

* initial guess
v.l(t) = vss;
pi_t.l(t) = pi_ss;
y.l(t) = yss;
z.l(t) = (1+rss)*(pi_t.l(t)/pi_ss)**phi_pi*(y.l(t)/yss)**phi_y - 1;
q.l(t) = ((1-theta*pi_t.l(t)**(alpha-1))/(1-theta))**(1/(1-alpha));
chi2.l(t) = chi2ss;
chi1.l(t) = q.l(t)*chi2.l(t)*(alpha-1)/alpha;
obj.l = 0;
epss.l(iii) = 0;

Option limrow=0,limcol=0,solprint=off;
OPTION ITERLIM = 500000;
OPTION RESLIM = 500000;
OPTION DOMLIM = 1000;
OPTION DNLP= conopt;
option decimals = 7;

model NewKeynesian /all/;

*******************************************************
* Step 2: Optimization step

set k grid space /1*11/;
alias(k,k2,d,d2);

* approximation domain
scalars
vmin    lower bound of v / 1 /
vmax    upper bound of v  / 1.045 /
bmin    lower bound of beta / 0.936 /
bmax    upper bound of beta / 1.056 / 
;

* Chebyshev parameters
parameters 
zv(k)
zb(k2)
extv
extb
extvmin
extbmin
extvmax
extbmax
dzv
dzb
;

zv(k) = - cos((2*ord(k)-1)/(2*card(k))*pi);
zb(k2) = - cos((2*ord(k2)-1)/(2*card(k2))*pi);
extv = -(zv('1')+1)*(vmax-vmin)/2/zv('1');
extb = -(zb('1')+1)*(bmax-bmin)/2/zb('1');
extvmin = vmin - extv;
extbmin = bmin - extb;
extvmax = vmax + extv;
extbmax = bmax + extb;
dzv = 2/(extvmax-extvmin);
dzb = 2/(extbmax-extbmin);

parameters
vinit(k)
binit(k2)
chi1sol(k,k2)
chi2sol(k,k2)
ysol(k,k2)
;

vinit(k) = extvmin + (1+zv(k))/dzv;
binit(k2) = extbmin + (1+zb(k2))/dzb;

loop((k,k2),
  v.fx('1') = vinit(k);

  betas('1') = binit(k2);
  loop(t$(ord(t) lt card(t)),
    betas(t+1) = exp((1-rho)*log(beta) + rho*log(betas(t)));
  );

  solve NewKeynesian using dnlp minimizing obj;

  chi1sol(k,k2) = chi1.l('1');
  chi2sol(k,k2) = chi2.l('1');
  ysol(k,k2) = y.l('1');
);


*******************************************************
* Step 3: Approximation step

parameters TTv(k,d)
        TTb(k2,d2)
;
TTv(k,'1') = 1;
TTv(k,'2') = zv(k);
loop(d$(ord(d)>=2 and ord(d)<card(d)),
  TTv(k,d+1) = 2*zv(k)*TTv(k,d) - TTv(k,d-1);
);
TTb(k2,'1') = 1;
TTb(k2,'2') = zb(k2);
loop(d2$(ord(d2)>=2 and ord(d2)<card(d2)),
  TTb(k2,d2+1) = 2*zb(k2)*TTb(k2,d2) - TTb(k2,d2-1);
);

parameters
coefsy(d,d2) 
coefschi1(d,d2)
coefschi2(d,d2)  
;

coefsy(d,d2) = 0;
coefsy(d,d2)$(ord(d)+ord(d2)<=card(d)+1) =
    4*sum((k,k2), TTv(k,d)*TTb(k2,d2)*ysol(k,k2)) / (card(k)*card(k2));
coefsy('1',d2) = coefsy('1',d2)/2;
coefsy(d,'1') = coefsy(d,'1')/2;

coefschi1(d,d2) = 0;
coefschi1(d,d2)$(ord(d)+ord(d2)<=card(d)+1) =
    4*sum((k,k2), TTv(k,d)*TTb(k2,d2)*chi1sol(k,k2)) / (card(k)*card(k2));
coefschi1('1',d2)= coefschi1('1',d2)/2;
coefschi1(d,'1') = coefschi1(d,'1')/2;

coefschi2(d,d2) = 0;
coefschi2(d,d2)$(ord(d)+ord(d2)<=card(d)+1) =
    4*sum((k,k2), TTv(k,d)*TTb(k2,d2)*chi2sol(k,k2)) / (card(k)*card(k2));
coefschi2('1',d2) = coefschi2('1',d2)/2;
coefschi2(d,'1') = coefschi2(d,'1')/2;

display coefsy, coefschi1, coefschi2;


*******************************************************
* Step 4: Error Checking

set k1 /1*10000/;
set k3 index of Gauss-Hermite quadrature rule /1*15/; 

Parameters
vmin1 	lower bound of v for error checking / 1 /
vmax1 	upper bound of v for error checking / 1.04 /
bmin1 	lower bound of beta for error checking / 0.96 /
bmax1 	upper bound of beta for error checking / 1.03 /
vv(k1) 
bb(k1) 
bbp(k1,k3)
vvp(k1)
ysol1(k1) 
ysol1p(k1,k3) 
chi1sol1(k1)
chi1sol1p(k1,k3)
chi2sol1(k1)
chi2sol1p(k1,k3)
qsol1(k1)
qsol1p(k1,k3)
pisol1(k1)
pisol1p(k1,k3)
zsol1(k1)
integrand1(k1,k3)
integrand2(k1,k3)
integrand3(k1,k3)
;


vv(k1) = Uniform(vmin1, vmax1);
bb(k1) = Uniform(bmin1, bmax1);

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


bbp(k1,k3) = exp((1-rho)*log(beta)+rho*log(bb(k1))+sigma*GHNodes(k3));

ysol1(k1) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefsy(d,d2)*cos((ord(d)-1)*arccos(dzv*(vv(k1)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(bb(k1)-extbmin)-1)) );
chi1sol1(k1) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi1(d,d2)*cos((ord(d)-1)*arccos(dzv*(vv(k1)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(bb(k1)-extbmin)-1)) );
chi2sol1(k1) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi2(d,d2)*cos((ord(d)-1)*arccos(dzv*(vv(k1)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(bb(k1)-extbmin)-1)) );

qsol1(k1) = chi1sol1(k1)*alpha/(alpha-1)/chi2sol1(k1);
pisol1(k1) = ((1 - qsol1(k1)**(1-alpha) * (1-theta))/theta )**(1/(alpha-1));
zsol1(k1) = (1+rss)*(pisol1(k1)/pi_ss)**phi_pi*(ysol1(k1)/yss)**phi_y - 1;
vvp(k1) = (1-theta)*qsol1(k1)**(-alpha)+theta*pisol1(k1)**alpha*vv(k1);

ysol1p(k1,k3) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefsy(d,d2)*cos((ord(d)-1)*arccos(dzv*(vvp(k1)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(bbp(k1,k3)-extbmin)-1)) );
chi1sol1p(k1,k3) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi1(d,d2)*cos((ord(d)-1)*arccos(dzv*(vvp(k1)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(bbp(k1,k3)-extbmin)-1)) );
chi2sol1p(k1,k3) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi2(d,d2)*cos((ord(d)-1)*arccos(dzv*(vvp(k1)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(bbp(k1,k3)-extbmin)-1)) );

qsol1p(k1,k3) = chi1sol1p(k1,k3)*alpha/(alpha-1)/chi2sol1p(k1,k3);
pisol1p(k1,k3) = ((1 - qsol1p(k1,k3)**(1-alpha) * (1-theta))/theta )**(1/(alpha-1));


integrand1(k1,k3) = bbp(k1,k3)*(1+max(0,zsol1(k1)))/pisol1p(k1,k3)*ysol1(k1)/ysol1p(k1,k3);
integrand2(k1,k3) = bbp(k1,k3)*pisol1p(k1,k3)**alpha*chi1sol1p(k1,k3);
integrand3(k1,k3) = bbp(k1,k3)*pisol1p(k1,k3)**(alpha-1)*chi2sol1p(k1,k3);


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

errs1(k1) = abs(sum(k3, GHWeights(k3)*integrand1(k1,k3)) - 1);
err1Linf = smax(k1, errs1(k1));
err1L1 = sum(k1, errs1(k1)) / card(k1);

errs2(k1) = abs(1 - (ysol1(k1)**(1+eta)*vvp(k1)**eta + theta*sum(k3, GHWeights(k3)*integrand2(k1,k3)))/chi1sol1(k1));
err2Linf = smax(k1, errs2(k1));
err2L1 = sum(k1, errs2(k1)) / card(k1);

errs3(k1) = abs(1 - (1/(1-sg) + theta*sum(k3, GHWeights(k3)*integrand3(k1,k3)))/chi2sol1(k1));
err3Linf = smax(k1, errs3(k1));
err3L1 = sum(k1, errs3(k1)) / card(k1);

display err1Linf, err1L1, err2Linf, err2L1, err3Linf, err3L1;

parameter numNeg number of negative z (so interest rate hits ZLB);

numNeg = 0;
loop(k1,
  if(zsol1(k1)<0,
    numNeg = numNeg+1;
  );
);

display numNeg;

*******************************************
* Report errors for plotting figures

set k4 /1*3/;
set k5 /1*100/;

Parameters
abb(k4)
avv(k5)  
abbp(k4,k3)
avvp(k4,k5)
aysol1(k4,k5) 
aysol1p(k4,k5,k3) 
achi1sol1(k4,k5)
achi1sol1p(k4,k5,k3)
achi2sol1(k4,k5)
achi2sol1p(k4,k5,k3)
aqsol1(k4,k5)
aqsol1p(k4,k5,k3)
apisol1(k4,k5)
apisol1p(k4,k5,k3)
azsol1(k4,k5)
arsol1(k4,k5)
aintegrand1(k4,k5,k3)
aintegrand2(k4,k5,k3)
aintegrand3(k4,k5,k3)
aerrs1(k4,k5)
aerrs2(k4,k5)
aerrs3(k4,k5)
aerrs(k4,k5)
;


avv(k5) = vmin1 + (ord(k5)-1)*(vmax1-vmin1)/(card(k5)-1);
abb('1') = bmin1;
abb('2') = beta;
abb('3') = bmax1;


abbp(k4,k3) = exp((1-rho)*log(beta)+rho*log(abb(k4))+sigma*GHNodes(k3));

aysol1(k4,k5) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefsy(d,d2)*cos((ord(d)-1)*arccos(dzv*(avv(k5)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(abb(k4)-extbmin)-1)) );
achi1sol1(k4,k5) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi1(d,d2)*cos((ord(d)-1)*arccos(dzv*(avv(k5)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(abb(k4)-extbmin)-1)) );
achi2sol1(k4,k5) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi2(d,d2)*cos((ord(d)-1)*arccos(dzv*(avv(k5)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(abb(k4)-extbmin)-1)) );

aqsol1(k4,k5) = achi1sol1(k4,k5)*alpha/(alpha-1)/achi2sol1(k4,k5);
apisol1(k4,k5) = ((1 - aqsol1(k4,k5)**(1-alpha) * (1-theta))/theta )**(1/(alpha-1));
azsol1(k4,k5) = (1+rss)*(apisol1(k4,k5)/pi_ss)**phi_pi*(aysol1(k4,k5)/yss)**phi_y - 1;
arsol1(k4,k5) = max(0, azsol1(k4,k5));
avvp(k4,k5) = (1-theta)*aqsol1(k4,k5)**(-alpha)+theta*apisol1(k4,k5)**alpha*avv(k5);

aysol1p(k4,k5,k3) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefsy(d,d2)*cos((ord(d)-1)*arccos(dzv*(avvp(k4,k5)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(abbp(k4,k3)-extbmin)-1)) );
achi1sol1p(k4,k5,k3) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi1(d,d2)*cos((ord(d)-1)*arccos(dzv*(avvp(k4,k5)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(abbp(k4,k3)-extbmin)-1)) );
achi2sol1p(k4,k5,k3) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi2(d,d2)*cos((ord(d)-1)*arccos(dzv*(avvp(k4,k5)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(abbp(k4,k3)-extbmin)-1)) );

aqsol1p(k4,k5,k3) = achi1sol1p(k4,k5,k3)*alpha/(alpha-1)/achi2sol1p(k4,k5,k3);
apisol1p(k4,k5,k3) = ((1 - aqsol1p(k4,k5,k3)**(1-alpha) * (1-theta))/theta )**(1/(alpha-1));


aintegrand1(k4,k5,k3) = abbp(k4,k3)*(1+arsol1(k4,k5))/apisol1p(k4,k5,k3)*aysol1(k4,k5)/aysol1p(k4,k5,k3);
aintegrand2(k4,k5,k3) = abbp(k4,k3)*apisol1p(k4,k5,k3)**alpha*achi1sol1p(k4,k5,k3);
aintegrand3(k4,k5,k3) = abbp(k4,k3)*apisol1p(k4,k5,k3)**(alpha-1)*achi2sol1p(k4,k5,k3);

aerrs1(k4,k5) = abs(sum(k3, GHWeights(k3)*aintegrand1(k4,k5,k3)) - 1);
aerrs2(k4,k5) = abs(1 - (aysol1(k4,k5)**(1+eta)*avvp(k4,k5)**eta + theta*sum(k3, GHWeights(k3)*aintegrand2(k4,k5,k3)))/achi1sol1(k4,k5));
aerrs3(k4,k5) = abs(1 - (1/(1-sg) + theta*sum(k3, GHWeights(k3)*aintegrand3(k4,k5,k3)))/achi2sol1(k4,k5));
aerrs(k4,k5) = max(aerrs1(k4,k5), max(aerrs2(k4,k5),aerrs3(k4,k5)));

file outputerr /newKeyn_errs.csv/;
outputerr.nw = 18;
outputerr.nr = 2;
outputerr.nz = 1e-15;
put outputerr;
outputerr.pc=5;
loop(k5,
  put avv(k5):14:6;
  loop(k4,
     put aerrs(k4,k5):14:6;
  );
  put /;
);


*******************************************
* Report solutions from impluse response anlysis for plotting figures

set k6 /1*21/;

Parameters
aabb(k6)
aavv(k6)
aaysol1(k6) 
aachi1sol1(k6)
aachi2sol1(k6)
aaqsol1(k6)
aapisol1(k6)
aazsol1(k6)
aarsol1(k6)
;


aavv('1') = vss;
aabb('1') = bmax1;
loop(k6$(ord(k6)<card(k6)),  
  aaysol1(k6) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefsy(d,d2)*cos((ord(d)-1)*arccos(dzv*(aavv(k6)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(aabb(k6)-extbmin)-1)) );
  aachi1sol1(k6) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi1(d,d2)*cos((ord(d)-1)*arccos(dzv*(aavv(k6)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(aabb(k6)-extbmin)-1)) );
  aachi2sol1(k6) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi2(d,d2)*cos((ord(d)-1)*arccos(dzv*(aavv(k6)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(aabb(k6)-extbmin)-1)) );

  aaqsol1(k6) = aachi1sol1(k6)*alpha/(alpha-1)/aachi2sol1(k6);
  aapisol1(k6) = ((1 - aaqsol1(k6)**(1-alpha) * (1-theta))/theta )**(1/(alpha-1));
  aazsol1(k6) = (1+rss)*(aapisol1(k6)/pi_ss)**phi_pi*(aaysol1(k6)/yss)**phi_y - 1;
  aarsol1(k6) = max(0, aazsol1(k6));

  aabb(k6+1) = exp((1-rho)*log(beta)+rho*log(aabb(k6)));
  aavv(k6+1) = (1-theta)*aaqsol1(k6)**(-alpha)+theta*aapisol1(k6)**alpha*aavv(k6);
);


file outputsol1 /newKeyn_sol1.csv/;
outputsol1.nw = 18;
outputsol1.nr = 2;
outputsol1.nz = 1e-15;
put outputsol1;
outputsol1.pc=5;
loop(k6$(ord(k6)<card(k6)),
  put aabb(k6):14:6;
  put aavv(k6):14:6;
  put aarsol1(k6):14:6;
  put aapisol1(k6):14:6;
  put aaysol1(k6):14:6;
  put /;
);


aavv('1') = vss;
aabb('1') = bmin1;
loop(k6$(ord(k6)<card(k6)),  
  aaysol1(k6) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefsy(d,d2)*cos((ord(d)-1)*arccos(dzv*(aavv(k6)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(aabb(k6)-extbmin)-1)) );
  aachi1sol1(k6) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi1(d,d2)*cos((ord(d)-1)*arccos(dzv*(aavv(k6)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(aabb(k6)-extbmin)-1)) );
  aachi2sol1(k6) = sum( (d,d2)$(ord(d)+ord(d2)<=card(d)+1),
        coefschi2(d,d2)*cos((ord(d)-1)*arccos(dzv*(aavv(k6)-extvmin)-1)) * cos((ord(d2)-1)*arccos(dzb*(aabb(k6)-extbmin)-1)) );

  aaqsol1(k6) = aachi1sol1(k6)*alpha/(alpha-1)/aachi2sol1(k6);
  aapisol1(k6) = ((1 - aaqsol1(k6)**(1-alpha) * (1-theta))/theta )**(1/(alpha-1));
  aazsol1(k6) = (1+rss)*(aapisol1(k6)/pi_ss)**phi_pi*(aaysol1(k6)/yss)**phi_y - 1;
  aarsol1(k6) = max(0, aazsol1(k6));

  aabb(k6+1) = exp((1-rho)*log(beta)+rho*log(aabb(k6)));
  aavv(k6+1) = (1-theta)*aaqsol1(k6)**(-alpha)+theta*aapisol1(k6)**alpha*aavv(k6);
);


file outputsol2 /newKeyn_sol2.csv/;
outputsol2.nw = 18;
outputsol2.nr = 2;
outputsol2.nz = 1e-15;
put outputsol2;
outputsol2.pc=5;
loop(k6$(ord(k6)<card(k6)),
  put aabb(k6):14:6;
  put aavv(k6):14:6;
  put aarsol1(k6):14:6;
  put aapisol1(k6):14:6;
  put aaysol1(k6):14:6;
  put /;
);

