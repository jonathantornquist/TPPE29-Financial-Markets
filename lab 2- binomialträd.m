%Löser alla uppgifter efter den implicita volatiliteten genom 
%att anropa ett underprogram kalla binomialtrad. 

%1 a) Bestäm den Implicita volatiliteten
r=-0.0005;
S_t=1915.220;
K=1900;
%Ska använda 3e fredagen i mars
T=(2+15*5-5)/252; %72
%Dagar i bruten vecka + veckor - helgdagar

%Löser ut den implicita volatiliteten med
%hjälp av intervall halvering
sigma=0;
sigma_ovre=1;

for n=1:50
d1=(log(S_t/K)+(r+(sigma^2/2))*T)/(sigma*sqrt(T));
d2=d1-(sigma*sqrt(T));
C=S_t*normcdf(d1)-K*exp(-r*T)*normcdf(d2);
 
if C>84.5
    sigma_ovre=sigma;
    sigma=(sigma_ovre+sigma_nedre)/2;
elseif C<84.5
    sigma_nedre=sigma;
    sigma=((sigma_ovre+sigma_nedre)/2);
else break 
end
end
sigma2=sigma;
fprintf("implicita volatiliteten =")
disp(sigma)
fprintf("Options pris via BSM, C=")
disp(C)


perioder=10;
utd=0;
optionstyp=0;
T_utd=0;
delta_T=T/perioder;
fprintf("Obligationspris via obligationsträd, 10 perioder=")
C1=binomialtrad(S_t,sigma,r,delta_T,K,perioder,utd,T_utd, optionstyp)
%Ska vara 84.5 enligt BSM och bild
fprintf("När perioder -> OO så går C->84.5, OK.    ");

% 1C) Plotta upp C för olika perioder (5-200)
Olika_C=zeros(1,195);
F=0;
utd=0;
optionstyp=0;
for m=5:200 
T=(2+15*5-5)/252; %72
T_utd=0;
perioder=m;
delta_T=T/perioder;
Olika_C(1,m-4)=binomialtrad(S_t,sigma,r,delta_T,K,perioder,utd,T_utd,optionstyp);

%Okej intervall [84.0775,84.9225]
if (84.0775<Olika_C(1,m-4) & F==0)
  fprintf("Värdet konvergerar +-0.5 procent efter")
   Antal_perioder=m
   F=1;
end
end
X=[5:200];
plot(X,Olika_C) %Olika C-värden för olika perioder. 


%Uppgift 2, givet
sigma=0.285;
utd=4.174;
K=250;
perioder=200;
T=(2+(41*5)-9)/252; %198
T_utd=(92-5)/252; %87
S_t=232.6;
%Vill ha 14.75
r=-0.0005;
S_t=S_t-utd*exp(-r*T_utd); %=S*
delta_T=T/perioder; 

fprintf("Estimerar C med rekombinerat binomialträd");
C=binomialtrad(S_t,sigma,r,delta_T,K,perioder,utd,T_utd,optionstyp)
fprintf("Avvikelse från faktiskt (C 14.75) ");
1-(C/14.75)

%Uppgift 3
sigma=sigma2;
utd=0;
optionstyp=1;
K=1900;
perioder=200;
T=(2+15*5-5)/252; %72
T_utd=0;
S_t=1915.22;
r=-0.0005;
delta_T=T/perioder;

fprintf("Estimerar C för down-and-out");
C=binomialtrad(S_t,sigma,r,delta_T,K,perioder,utd,T_utd,optionstyp)


%Vårt underprogram som vi använder för binomialträd
function C_i= binomialtrad(S_t,sigma,r,delta_T,K,perioder,utd,T_utd,optionstyp)

u=exp(sigma*sqrt(delta_T));
d=exp(-sigma*sqrt(delta_T));
q=(exp(r*delta_T)-d)/(u-d);
aktietrad=zeros(perioder,perioder);
aktietrad(1,1)=S_t;


for k=2:perioder;
    aktietrad(k,k)=aktietrad(k-1,k-1)*d; %Diagonalen, bara nedgångar. 
end

raknare=2; %Håller koll på kolumn, ska flyttas oss en till höger.
for i=1:perioder %Bara uppgångar
    for j=raknare:perioder %(rad,kolumn)
        aktietrad(i,j)=aktietrad(i,j-1)*u;    
end
raknare=raknare+1; 
end
aktietrad;

optionstrad=zeros(perioder,perioder);
for u=1:perioder  %Löser sista kolumnen
    if optionstyp==0 %köp
    optionstrad(u,perioder)=aktietrad(u,perioder)-K;
    elseif optionstyp==1 %sälj och down-and-out
    optionstrad(u,perioder)=K-aktietrad(u,perioder);
    if aktietrad(u,perioder)<S_t*0.92 %Om aktiepriset faller mer än 8% nollas optionsvärdet i noden
       optionstrad(u,perioder)=0;
    end
    end
if optionstrad(u,perioder)<0
   optionstrad(u,perioder)=0;
end
end


%Rekursiv lösning av binomialträd.
k=1;
for j=perioder-1:-1:1 %Loop i reverse
utd2=utd*exp(-r*delta_T*j); %nuvärdet av utdelningen
for i=1:perioder-k %Löser det för kolumn j

if utd>0  & j*delta_T<=T_utd %tar max(C_i,j;S_i,j+DIV-K),   
    optionstrad(i,j)=max(exp(-r*delta_T)*(q*optionstrad(i,j+1)+(1-q)*optionstrad(i+1,j+1)),aktietrad(i,j)+utd2-K); 
elseif  utd>0 & j*delta_T>T_utd %Om vi får utdelning eller inte. Om passerat utd datum. 
optionstrad(i,j)=max(exp(-r*delta_T)*(q*optionstrad(i,j+1)+(1-q)*optionstrad(i+1,j+1)),aktietrad(i,j)-K); 
elseif optionstyp==1
    optionstrad(i,j)=exp(-r*delta_T)*(q*optionstrad(i,j+1)+(1-q)*optionstrad(i+1,j+1));
    if aktietrad(i,j)<S_t*0.92
       optionstrad(i,j)=0;
    end
else
   %I fallen när utdelningen är noll. (Dvs uppg.1)
    optionstrad(i,j)=exp(-r*delta_T)*(q*optionstrad(i,j+1)+(1-q)*optionstrad(i+1,j+1));
end


 
end
k=k+1;
end

C_i=optionstrad(1,1);
optionstrad;
end 
