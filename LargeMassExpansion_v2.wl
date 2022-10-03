(* ::Package:: *)

<<LiteRed`; (* Loading package *)
Remove[Global`IBPReduce]
Remove[Global`j]
Remove[Global`Dinv]
SetDirectory[NotebookDirectory[]];
ClearAll["d",s,t,"u","mm","a","y"]
SetDim[d]; (* Set dimensionality*)
Declare[{k,q1,q2,q3},Vector,{s,t,mm},Number];
sp[q1]=sp[q2]=sp[q3]=0; (* constrains *)
sp[q1+q2]=s;
sp[q1+q3]=t;
sp[q2+q3]=-t-s;
sp[q1,q2]=s/2;
sp[q1,q3]=t/2;
sp[q2,q3]=(-t-s)/2;

NewBasis[Base,{mm-sp[k],mm-sp[k+q1],mm-sp[k+q1+q2],-sp[k+q1+q2+q3]},{k},Directory->"Base"]; (* Set up base*)
GenerateIBP[Base]; (*Generate IBP Relations*)
AnalyzeSectors[Base];
FindSymmetries[Base];
SolvejSector/@UniqueSectors[Base];

nmaximal=10; (* Set max order in (s/mm) *)
epsmaximal=1; (* Set max order in epsilon *) (*min 1*)

MasInt={j[Base,0,0,1,0]->mm*scale*Exp[EulerGamma*eps]*Gamma[1-d/2]}; (*Known Result for Tadpole*)

SolvejDGL2[j_,x_,nmaxi_,epsmaxi_,bc__]:=Module[{yy,coeffs,DGL,CoeffList1,CoeffList2,SeriesDGL1,SeriesDGL2}, (*Input: MasInt,kinematic variable (e.g "s"), max orders, List of boundaries *)
yy=Sum[eps^m*(Sum[a[m,n]*(x/mm)^n,{n,0,nmaxi}]+Warn[f,j]*x^(nmaxi+1)),{m,-1,epsmaxi}]+Warn[e]*eps^(epsmaxi+1); (* Ansatz *)
coeffs={}; (* Set up empty coefficient List *)
DGL=D[yy,{x,1}]-IBPReduce[Dinv[j,x]]/.MasInt /.{j->yy,d->4-2eps}; (* Get DEq and insert known Integrals *)
SeriesDGL1=Normal[Series[DGL,{x,0,nmaxi-1}]];
SeriesDGL2=Normal[Series[SeriesDGL1,{eps,0,epsmaxi}]]; (* Series Expansion in s and epsilon *)
CoeffList1=CoefficientList[eps*SeriesDGL2,eps];
CoeffList2=CoefficientList[x*CoeffList1,x];  (* Build up Coefficient Lists *)
For[i=1,i<Length[CoeffList2]+1,i++,
AppendTo[coeffs,First[Solve[CoeffList2[[i]]==0,MapThread[a,{ConstantArray[-2+i,nmaxi+1],Range[0,nmaxi]}]]]];
]; (* Comparision of coefficients, order by order *)
For[i=2,i<Length[coeffs]+1,i++,
coeffs[[i]]=coeffs[[i]]/.Flatten[coeffs[[1;;i-1]]] 
]; (* Insert solutions step by step *)
coeffs=Flatten[coeffs];
yy/.coeffs/.bc (* Output: Ansatz with calculated coefficients and boundary conditions*)
]

AppendToMasInt[j_,s_,nmaxi_,epsmaxi_,bc__]:=Module[{yysol}, (* Append results to list of masters*)
yysol=SolvejDGL2[j,s,nmaxi,epsmaxi,bc];
If[FilterRules[MasInt,j]=={},AppendTo[MasInt,j->yysol],Print["F\[UDoubleDot]r das Masterintegral ",j," ist schon eine Regel vorhanden."]]

]
f[n_,m_]:=n->m
listToRules[l1__,l2__]:=MapThread[f,{l1,l2}]

BC01[nmax_]:=Module[{temp}, (*Calculate boundary conditons from Two-Point-Integral with s\[Rule]0 *)
temp=CoefficientList[eps*Series[scale*Exp[EulerGamma*eps]*Gamma[eps],{eps,0,nmax}]*Assuming[t<0,Integrate[Series[((1-x-y*x*(x-1)))^(-eps),{eps,0,nmax+1}],{x,0,1}]],eps];
listToRules[a[#,0]&/@Range[-1,nmax],temp/.y->t/mm]
]
BC1[nmax_]:=Module[{temp}, (*Calculate boundary conditons from Three-Point-Integral with s\[Rule]0 *)
temp=(mm)^(-1)*Normal[Integrate[Normal[Integrate[1/B*Series[scale*Exp[EulerGamma*eps]*Gamma[1+eps]*w^(-1-eps),{eps,0,nmax}],{w,A,A+B}]]/.A->(1-z)/.B->(x*z),{z,0,1}]];
listToRules[a[#,0]&/@Range[-1,nmax],CoefficientList[eps*temp,eps]/.x->t/mm]
]

Timing[AppendToMasInt[j[Base,0,1,0,1],s,nmaximal,epsmaximal,BC01[epsmaximal]];
AppendToMasInt[j[Base,0,1,1,1],s,nmaximal,epsmaximal,BC1[epsmaximal]];
AppendToMasInt[j[Base,1,1,0,1],s,nmaximal,epsmaximal,BC1[epsmaximal]];
AppendToMasInt[j[Base,1,0,1,0],s,nmaximal,epsmaximal,{}];
AppendToMasInt[j[Base,1,0,1,1],s,nmaximal,epsmaximal,{}];
AppendToMasInt[j[Base,1,1,1,0],s,nmaximal,epsmaximal,{}];
AppendToMasInt[j[Base,1,1,1,1],s,nmaximal,epsmaximal,{}];
(*MasInt=MasInt/.scale->1/.mm/mu2->1;*)
Export["MasIntScale2.m",MasInt]]
















