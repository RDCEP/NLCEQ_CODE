*----------------------------------------------------------------------
* This program solves the 2-country RBC model using NLCEQ
*
* Authors: Yongyang Cai, University of Chicago and Hoover Institution
*	   Jevgenijs Steinbuks, The World Bank
*
* If using material from this code, the user should cite the following paper:
* Cai, Y., K.L. Judd, and J. Steinbuks (2016). A nonlinear certainty equivalent* approximation method for dynamic stochastic problems. Quantitative Economics (forthcoming).
*----------------------------------------------------------------------

scalars
beta0 discount rate /0.99/
alpha capital cost share /0.36/
delta capital stock depreciation /0.025/
rho persistence parameter for TFP growth /0.95/
phi adjustment cost parameter /0.5/
sigma standard deviation of the stochastic productivity level / 0.01 /
gamma intertemporal elasticity of substitution /0.5/
gammahat utility parameter
eta Frisch elasticity of labor supply /0.5/
etahat utility parameter
A technology parameter
B relative weight of consumption and leisure
;

A = (1 - beta0) / (alpha * beta0);
gammahat = 1-(1/gamma);
B = (1 - alpha) * (A ** gammahat);
etahat = 1+(1/eta);


set j agents /1*2/;
alias(j,j1);

parameter tau(j) Negishi weight;
tau(j) = 1/card(j);


*************************************************
* Step 1: Define the transformed deterministic problem 

set t time /1*51/;

Parameters
beta(t) discount rate
theta(j,t) productivity shock
;

beta(t) = beta0**(ord(t) - 1);
theta(j,'1') = 1;
loop (t$(ord(t) lt card(t)),
  theta(j,t+1) = theta(j,t)**rho;
);


Variables
obj objective criterion
u(j,t) utility
i(j,t) investment
;

Positive variables
k(j,t) capital stock
c(j,t) consumption
l(j,t) labor supply
y(j,t) output
;

Equations
obj1 Objective function
util(j,t) Utility function
output(j,t) Output
cap(j,t) Law of Motion for Capital Stock
bc(t) budget constraint
;

cap(j,t)$(ord(t) lt card(t)) .. 
k(j,t+1) =e= (1-delta)*k(j,t) + i(j,t);

output(j,t)$(ord(t) lt card(t)) .. 
y(j,t) =e= theta(j,t)*A*(k(j,t)**alpha) * (l(j,t)**(1-alpha));

bc(t)$(ord(t) lt card(t)) .. 
sum(j, c(j,t) + i(j,t)) =e= sum(j, y(j,t) + delta * k(j,t)) - sum(j, (phi/2)*k(j,t)*sqr(i(j,t)/k(j,t)-delta));

util(j,t)$(ord(t) lt card(t)) .. 
u(j,t) =e= (c(j,t)**gammahat)/gammahat - B * (l(j,t)**etahat)/etahat;

obj1 .. 
obj =e= sum(j, tau(j) * sum(t$(ord(t) lt card(t)), beta(t)*u(j,t))) + 
  sum(j, tau(j) * sum(t$(ord(t)=card(t)), beta(t)*((( (A*(k(j,t)**alpha))**gammahat )/gammahat-B)/(1-beta0)))) ;


* Bound Constraints

k.lo(j,t) = 0.01;
k.up(j,t) = 100;
c.lo(j,t) = 0.001;
c.up(j,t) = 100;
l.lo(j,t) = 0.001;
l.up(j,t) = 100;

* Initial Guess

i.l(j,t) = delta;
k.l(j,t) = 1;
l.l(j,t) = 1;
y.l(j,t) = theta(j,t)*A*(k.l(j,t)**alpha)*(l.l(j,t)**(1-alpha));
c.l(j,t) = A;
u.l(j,t) =(c.l(j,t)**gammahat)/gammahat - B * (l.l(j,t)**etahat)/etahat;
obj.l = sum(j, tau(j) * sum(t$(ord(t) lt card(t)), beta(t)*u.l(j,t))) +
  sum(j, tau(j) * sum(t$(ord(t)=card(t)), beta(t)*((((A*(k.l(j,t)**alpha))**gammahat)/gammahat-B)/(1-beta0))));


Option limrow=0,limcol=0,solprint=off;
OPTION ITERLIM = 500000;
OPTION RESLIM = 500000;
OPTION DOMLIM = 1000;
OPTION NLP= conopt;
option decimals = 7;

model busc /all/;


*******************************************************
* Step 2: Optimization step

set v1 grid space /1*5/;
alias(v1,v2,v3,v4,d1,d2,d3,d4);

* approximation domain
scalars
kmin    lower bound of k  / 0.5 /
kmax    upper bound of k  / 1.5 /
amin    lower bound of theta  / 0.5 /
amax    upper bound of theta  / 1.5 /
;

* Chebyshev parameters

parameters zk(v1)
        za(v3)
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
za(v3) = - cos((2*ord(v3)-1)/(2*card(v3))*pi);
extk = -(zk('1')+1)*(kmax-kmin)/2/zk('1');
exta = -(za('1')+1)*(amax-amin)/2/za('1');
extkmin = kmin - extk;
extamin = amin - exta;
extkmax = kmax + extk;
extamax = amax + exta;
dzk = 2/(extkmax-extkmin);
dza = 2/(extamax-extamin);

Parameters
kinit(v1) Grid: Capital Stock of Country
ainit(v3) Grid: Productivity of Country 
;

kinit(v1) = extkmin + (1+zk(v1))/dzk;
ainit(v3) = extamin + (1+za(v3))/dza;

Parameters
solC(j,v1,v2,v3,v4) consumption
solI(j,v1,v2,v3,v4) investment
solL(j,v1,v2,v3,v4) labor supply
;

loop ((v1,v2,v3,v4),

  k.fx('1','1') = kinit(v1);
  k.fx('2','1') = kinit(v2);
  theta('1','1') = ainit(v3);
  theta('2','1') = ainit(v4);

  loop (t$(ord(t) lt card(t)),
    theta(j,t+1) = theta(j,t)**rho;
  );

  solve busc using nlp maximizing obj;

  solC(j,v1,v2,v3,v4) = c.l(j,'1');
  solI(j,v1,v2,v3,v4) = i.l(j,'1');
  solL(j,v1,v2,v3,v4) = l.l(j,'1');
);


************************************************
* Step 3: Approximation step

parameters TTk(v1,d1)
        TTa(v3,d3)
;
TTk(v1,'1') = 1;
TTk(v1,'2') = zk(v1);
loop(d1$(ord(d1)>=2 and ord(d1)<card(d1)),
  TTk(v1,d1+1) = 2*zk(v1)*TTk(v1,d1) - TTk(v1,d1-1);
);
TTa(v3,'1') = 1;
TTa(v3,'2') = za(v3);
loop(d3$(ord(d3)>=2 and ord(d3)<card(d3)),
  TTa(v3,d3+1) = 2*za(v3)*TTa(v3,d3) - TTa(v3,d3-1);
);

parameters 
coefsC(j,d1,d2,d3,d4)
coefsI(j,d1,d2,d3,d4)
coefsL(j,d1,d2,d3,d4)
;

coefsC(j,d1,d2,d3,d4) = 0;
coefsC(j,d1,d2,d3,d4)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)<=card(d1)+3) =
    16*sum((v1,v2,v3,v4), TTk(v1,d1)*TTk(v2,d2)*TTa(v3,d3)*TTa(v4,d4)*solC(j,v1,v2,v3,v4)) / (card(v1)*card(v2)*card(v3)*card(v4));
coefsC(j,'1',d2,d3,d4) = coefsC(j,'1',d2,d3,d4)/2;
coefsC(j,d1,'1',d3,d4) = coefsC(j,d1,'1',d3,d4)/2;
coefsC(j,d1,d2,'1',d4) = coefsC(j,d1,d2,'1',d4)/2;
coefsC(j,d1,d2,d3,'1') = coefsC(j,d1,d2,d3,'1')/2;

coefsI(j,d1,d2,d3,d4) = 0;
coefsI(j,d1,d2,d3,d4)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)<=card(d1)+3) =
    16*sum((v1,v2,v3,v4), TTk(v1,d1)*TTk(v2,d2)*TTa(v3,d3)*TTa(v4,d4)*solI(j,v1,v2,v3,v4)) / (card(v1)*card(v2)*card(v3)*card(v4));
coefsI(j,'1',d2,d3,d4) = coefsI(j,'1',d2,d3,d4)/2;
coefsI(j,d1,'1',d3,d4) = coefsI(j,d1,'1',d3,d4)/2;
coefsI(j,d1,d2,'1',d4) = coefsI(j,d1,d2,'1',d4)/2;
coefsI(j,d1,d2,d3,'1') = coefsI(j,d1,d2,d3,'1')/2;

coefsL(j,d1,d2,d3,d4) = 0;
coefsL(j,d1,d2,d3,d4)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)<=card(d1)+3) =
    16*sum((v1,v2,v3,v4), TTk(v1,d1)*TTk(v2,d2)*TTa(v3,d3)*TTa(v4,d4)*solL(j,v1,v2,v3,v4)) / (card(v1)*card(v2)*card(v3)*card(v4));
coefsL(j,'1',d2,d3,d4) = coefsL(j,'1',d2,d3,d4)/2;
coefsL(j,d1,'1',d3,d4) = coefsL(j,d1,'1',d3,d4)/2;
coefsL(j,d1,d2,'1',d4) = coefsL(j,d1,d2,'1',d4)/2;
coefsL(j,d1,d2,d3,'1') = coefsL(j,d1,d2,d3,'1')/2;

display coefsC, coefsI, coefsL;



*******************************************************
* Step 4: Error checking

set k1 /1*10000/;

set ii indices of Gaussian Hermite quadrature nodes /1*7/;
alias (ii, ii2);

parameter GHNodes(ii) Gaussian Hermite quadrature nodes;
GHNodes('1') =  -2.651961356835233*sqrt(2.0);
GHNodes('2') =  -1.673551628767471*sqrt(2.0);
GHNodes('3') =  -0.8162878828589647*sqrt(2.0);
GHNodes('4') =  0;
GHNodes('5') =  0.8162878828589647*sqrt(2.0);
GHNodes('6') =  1.673551628767471*sqrt(2.0);
GHNodes('7') =  2.651961356835233*sqrt(2.0);

parameter GHWeights(ii) Gaussian Hermite quadrature weights;
GHWeights('1') =  0.0009717812450995192/sqrt(pi);
GHWeights('2') =  0.05451558281912703/sqrt(pi);
GHWeights('3') =  0.4256072526101278/sqrt(pi);
GHWeights('4') =  0.8102646175568073/sqrt(pi);
GHWeights('5') =  0.4256072526101278/sqrt(pi);
GHWeights('6') =  0.05451558281912703/sqrt(pi);
GHWeights('7') =  0.0009717812450995192/sqrt(pi);

parameter GHWeights2D(ii,ii2) Gaussian Hermite quadrature weights in 2D;
GHWeights2D(ii,ii2) = GHWeights(ii)* GHWeights(ii2);

parameter chol(j,j1) cholesky factorization of correlation matrix of (e1+e) and (e2+e);
chol(j,j1) = 0;
chol('1','1') = 1;
chol('2','1') = 0.5;
chol('2','2') = 0.866;

* error checking domain
scalars
kmin1 minimal capital / 0.7 /
kmax1 maximal capital / 1.3 /
amin1 minimal productivity / 0.7 /
amax1 maximal productivity / 1.3 /
;


Parameters
kk(j,k1) capital stock
aa(j,k1) Stochastic productivity
aap1(k1,ii) Stochastic realization of productivity shock
aap2(k1,ii,ii2) Stochastic realization of productivity shock
kpsol1(j,k1)
csol1(j,k1)
csol1p(j,k1,ii,ii2)
Isol1(j,k1)
Isol1p(j,k1,ii,ii2)
Lsol1(j,k1)
Lsol1p(j,k1,ii,ii2)
integrand(j,k1,ii,ii2)
;

kk(j,k1) = Uniform(kmin1, kmax1);
aa(j,k1) = Uniform(amin1, amax1);

* ln(theta_j^+) = rho*ln(theta_j) + sigma * (eps_j + eps)
aap1(k1,ii) = aa('1',k1)**rho * exp(sigma*sqrt(2)*GHNodes(ii));
aap2(k1,ii,ii2) = aa('2',k1)**rho * exp(sigma*sqrt(2)*(chol('2','1')*GHNodes(ii)+chol('2','2')*GHNodes(ii2)));

csol1(j,k1) = sum( (d1,d2,d3,d4)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)<=card(d1)+3), coefsC(j,d1,d2,d3,d4)*
	cos((ord(d1)-1)*arccos(dzk*(kk('1',k1)-extkmin)-1)) * 
	cos((ord(d2)-1)*arccos(dzk*(kk('2',k1)-extkmin)-1)) * 
	cos((ord(d3)-1)*arccos(dza*(aa('1',k1)-extamin)-1)) * 
      	cos((ord(d4)-1)*arccos(dza*(aa('2',k1)-extamin)-1)) ); 

Isol1(j,k1) = sum( (d1,d2,d3,d4)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)<=card(d1)+3), coefsI(j,d1,d2,d3,d4)*
	cos((ord(d1)-1)*arccos(dzk*(kk('1',k1)-extkmin)-1)) * 
	cos((ord(d2)-1)*arccos(dzk*(kk('2',k1)-extkmin)-1)) * 
	cos((ord(d3)-1)*arccos(dza*(aa('1',k1)-extamin)-1)) * 
      	cos((ord(d4)-1)*arccos(dza*(aa('2',k1)-extamin)-1)) ); 

Lsol1(j,k1) = sum( (d1,d2,d3,d4)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)<=card(d1)+3), coefsL(j,d1,d2,d3,d4)*
	cos((ord(d1)-1)*arccos(dzk*(kk('1',k1)-extkmin)-1)) * 
	cos((ord(d2)-1)*arccos(dzk*(kk('2',k1)-extkmin)-1)) * 
	cos((ord(d3)-1)*arccos(dza*(aa('1',k1)-extamin)-1)) * 
      	cos((ord(d4)-1)*arccos(dza*(aa('2',k1)-extamin)-1)) ); 

kpsol1(j,k1) = (1-delta)*kk(j,k1) + Isol1(j,k1);

csol1p(j,k1,ii,ii2) = sum( (d1,d2,d3,d4)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)<=card(d1)+3), coefsC(j,d1,d2,d3,d4)*
	cos((ord(d1)-1)*arccos(dzk*(kpsol1('1',k1)-extkmin)-1)) * 
	cos((ord(d2)-1)*arccos(dzk*(kpsol1('2',k1)-extkmin)-1)) * 
	cos((ord(d3)-1)*arccos(dza*(aap1(k1,ii)-extamin)-1)) * 
      	cos((ord(d4)-1)*arccos(dza*(aap2(k1,ii,ii2)-extamin)-1)) ); 

Isol1p(j,k1,ii,ii2) = sum( (d1,d2,d3,d4)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)<=card(d1)+3), coefsI(j,d1,d2,d3,d4)*
	cos((ord(d1)-1)*arccos(dzk*(kpsol1('1',k1)-extkmin)-1)) * 
	cos((ord(d2)-1)*arccos(dzk*(kpsol1('2',k1)-extkmin)-1)) * 
	cos((ord(d3)-1)*arccos(dza*(aap1(k1,ii)-extamin)-1)) * 
      	cos((ord(d4)-1)*arccos(dza*(aap2(k1,ii,ii2)-extamin)-1)) ); 

Lsol1p(j,k1,ii,ii2) = sum( (d1,d2,d3,d4)$(ord(d1)+ord(d2)+ord(d3)+ord(d4)<=card(d1)+3), coefsL(j,d1,d2,d3,d4)*
	cos((ord(d1)-1)*arccos(dzk*(kpsol1('1',k1)-extkmin)-1)) * 
	cos((ord(d2)-1)*arccos(dzk*(kpsol1('2',k1)-extkmin)-1)) * 
	cos((ord(d3)-1)*arccos(dza*(aap1(k1,ii)-extamin)-1)) * 
      	cos((ord(d4)-1)*arccos(dza*(aap2(k1,ii,ii2)-extamin)-1)) ); 

integrand('1',k1,ii,ii2) = csol1p('1',k1,ii,ii2)**(-1/gamma) * 
	( 1+(phi/2)*(Isol1p('1',k1,ii,ii2)/kpsol1('1',k1)-delta) * (2-delta+Isol1p('1',k1,ii,ii2)/kpsol1('1',k1)) + 
	  aap1(k1,ii)*A*alpha * kpsol1('1',k1)**(alpha-1) * Lsol1p('1',k1,ii,ii2)**(1-alpha) );

integrand('2',k1,ii,ii2) = csol1p('2',k1,ii,ii2)**(-1/gamma) * 
	( 1+(phi/2)*(Isol1p('2',k1,ii,ii2)/kpsol1('2',k1)-delta) * (2-delta+Isol1p('2',k1,ii,ii2)/kpsol1('2',k1)) + 
	  aap2(k1,ii,ii2)*A*alpha * kpsol1('2',k1)**(alpha-1) * Lsol1p('2',k1,ii,ii2)**(1-alpha) );


parameters
errs1(j,k1)
errs2(k1)
errs3(j,k1)
errs4(k1)
err1Linf
err2Linf
err3Linf
err4Linf
err1L1
err2L1
err3L1
err4L1
;

errs1(j,k1) = abs( ( beta0*sum((ii,ii2), GHWeights2D(ii,ii2)*integrand(j,k1,ii,ii2)) ) / ( csol1(j,k1)**(-1/gamma) * (1+phi*(Isol1(j,k1)/kk(j,k1)-delta)) ) - 1 );
err1Linf = smax((j,k1), errs1(j,k1));
err1L1 = sum((j,k1), errs1(j,k1)) / (card(j)*card(k1));

errs2(k1) = abs( csol1('2',k1)**(-1/gamma) / csol1('1',k1)**(-1/gamma) * (tau('2')/tau('1')) - 1 );
err2Linf = smax(k1, errs2(k1));
err2L1 = sum(k1, errs2(k1)) / card(k1);

errs3(j,k1) = abs( csol1(j,k1)**(-1/gamma) * aa(j,k1) * (A*(1-alpha) * kk(j,k1)**alpha * Lsol1(j,k1)**(-alpha)) / (-B*Lsol1(j,k1)**(1/eta)) + 1 );
err3Linf = smax((j,k1), errs3(j,k1));
err3L1 = sum((j,k1), errs3(j,k1)) / (card(j)*card(k1));

errs4(k1) = abs( sum(j, csol1(j,k1)+Isol1(j,k1)-delta*kk(j,k1)+(phi/2)*kk(j,k1)*sqr(Isol1(j,k1)/kk(j,k1)-delta)) / sum(j,aa(j,k1)*A*kk(j,k1)**alpha*Lsol1(j,k1)**(1-alpha)) - 1 );
err4Linf = smax(k1, errs4(k1));
err4L1 = sum(k1, errs4(k1)) / card(k1);

display err1Linf, err1L1, err2Linf, err2L1;
display err3Linf, err3L1, err4Linf, err4L1;
