(* ::Package:: *)

<<LiteRed`; (* Loading package *)
Remove[Global`IBPReduce]
Remove[Global`j]
Remove[Global`Dinv]
SetDirectory[NotebookDirectory[]];
ClearAll["d",s,t,"u","mm","a","y"]
SetDim[d];
Declare[{k,q1,q2,q3},Vector,{s,t,mm},Number];
sp[q1]=sp[q2]=sp[q3]=0; (* constrains *)
sp[q1+q2]=s;
sp[q1+q3]=t;
sp[q2+q3]=-t-s;
sp[q1,q2]=s/2;
sp[q1,q3]=t/2;
sp[q2,q3]=(-t-s)/2;
NewBasis[Base,{mm-sp[k],mm-sp[k+q1],mm-sp[k+q1+q2],-sp[k+q1+q2+q3]},{k},Directory->"Base"]; (* Set up base*)
GenerateIBP[Base];
AnalyzeSectors[Base];
FindSymmetries[Base];
SolvejSector/@UniqueSectors[Base];
nmaximal=10; (* set up max order of (mm/s) or (mm/u)*)
(* eps order fixed due to limited order in boundaries*)
MasInt2={j[Base,0,0,1,0]->mm*(mm/mu2)^(-eps)*Exp[EulerGamma*eps]*Gamma[-1+eps]}; (* Known Tadpole *)

FF[nn_,nmax_]:=Module[{temp=0,yy0}, (* Ansatz *)
yy0=1/eps*Sum[(A[n]+H[n]*Log[-x])*x^n,{n,0,nmax}]+Warn[0]*x^(nmax+1)/eps;
For[i=1,i<nn+1,i++,
yy0=yy0+eps^(i-1)*Sum[(B[i-1,n]+Sum[C[i-1,n,m]*Log[-x]^m,{m,1,i+1}])*x^n,{n,0,nmax}]+Warn[i]*x^(nmax+1)
];
yy0+Warn[epsi]*eps^(nn)
]

SolvejDGLm3[j_,var_,nmax_,epsmax_,bc_]:=Module[{yy,yy1,yy2,yysol,coeffs,DGL,CoeffList1,CoeffList2,CoeffList3,SeriesDGL1,SeriesDGL2}, (* Similar routine to LargeMassExpansion*)
yy2=FF[epsmax+1,nmax]/.x->(mm/var);
DGL=D[yy2,mm]-IBPReduce[Dinv[j,mm]]/.MasInt2/.{j->yy2,d->4-2eps,t->-s-u};
SeriesDGL1=Normal[Series[DGL,{mm,0,nmax-1}]];
SeriesDGL2=Normal[Series[SeriesDGL1,{eps,0,epsmax}]];
CoeffList1=CoefficientList[eps*SeriesDGL2,eps];
CoeffList2=CoefficientList[mm*CoeffList1,mm]/.(Log[mm/mu2])->(Log[mm]-Log[mu2])/.Log[-mm/var]->(Log[mm]-Log[-var]);
CoeffList3=Flatten[CoefficientList[CoeffList2,Log[mm]],1]/.var->T[w];
If[StringMatchQ[ToString[j],"j[Base, 1, 1, 1, 1]"],coeffs=First[Solve[CoeffList3[[1]]==0]];
For[i=2,i<Length[CoeffList3]+1,i++,
yy1=CoeffList3[[i]]/.coeffs;
AppendTo[coeffs,Flatten[Solve[yy1==0],1]];
coeffs=Flatten[coeffs];
],coeffs=Flatten[Solve[CoeffList3==0],1]/.T[w]->var];
yysol=Collect[yy2/.coeffs/.bc,x]/.T[w]->var/.x->mm/var/.t->-s-u/.(Log[mm/mu2])->(Log[mm]-Log[mu2])/.Log[-mm/var]->(Log[mm]-Log[-var])]

Boundaries=Import["/Users/lennard/Desktop/LoopCalc/boundaries-small_m.m"]/.ep->eps/.m^2->mm/.m2->mm; (* Extract Data and modify it*)

GetBoundaries2[j_,var_,epsmax_,BC_]:=Module[{ListA,ListB,coeffs},
coeffs=Flatten[Solve[CoefficientList[eps*SolvejDGLm3[j,var,0,epsmax,{}]/.Warn[n_]->0,{eps,Log[-mm/var]}]==CoefficientList[eps*j/.BC,{eps,Log[-mm/var]}]]] (* Compare Boundaries with routine result for s,u\[Rule]0 *)]

GetBoundaries3[j_,var_,epsmax_,bc_]:=GetBoundaries3[j,var,epsmax,bc]=Module[{coeffs,ListA},
ListA=CoefficientList[eps*SolvejDGLm3[j,var,0,epsmax,{}]/.Warn[n_]->0,{eps,Log[mm]}]-CoefficientList[eps*j/.bc/.(Log[mm/mu2])->(Log[mm]-Log[mu2])/.Log[-mm/var]->(Log[mm]-Log[-var]),{eps,Log[mm]}];
ListA=ListA/.Log[mu2]->0/.Log[-var]->Log[-var/mu2]//FunctionExpand;
coeffs=Flatten[Solve[ListA==0]]
]

GetBoundaries4[bc_]:=GetBoundaries4[bc]=Module[{coeffs,ListA,ListB},
ListA=CoefficientList[eps*SolvejDGLm3[j[Base,1,1,1,1],u,0,0,{}]/.Warn[n_]->0,{eps,Log[mm]}]/.Log[-s/mu2]->(Log[-s]-Log[mu2])//Expand;
ListB=CoefficientList[eps*j[Base,1,1,1,1]/.Boundaries/.(Log[mm/mu2])->(Log[mm]-Log[mu2])/.Log[-mm/u]->(Log[mm]-Log[-u]),{eps,Log[mm]}];
coeffs=Flatten[Normal[Solve[ListA==ListB]]]/.Log[s/u]->(Log[-s]-Log[-u])/.Log[-u]->0/.Log[-s]->Log[s/u]//Expand]

Rule2s={Log[mm]-Log[-s]->Log[-mm/s],Log[-s]->Log[-s/mu2]+Log[mu2]};
Rule2u={Log[mm]->Log[-mm/u]+Log[-u],Log[mu2]->Log[-mu2/u]+Log[-u],Log[-mu2/u]->-Log[-u/mu2]};
Rule3s={Log[mm]-Log[-s]->Log[-mm/s],Log[mu2]->Log[-mu2/s]+Log[-s],Log[-s/mu2]->-Log[-mu2/s],Log[-s/mu2]->-Log[-mu2/s]};
Rule3u={Log[mm]-Log[-u]->Log[-mm/u],Log[mu2]->Log[-mu2/u]+Log[-u],Log[-mu2/u]->-Log[-u/mu2]};
RuleBox={Log[mm]-Log[-u]->Log[-mm/u],Log[-u/mu2]->Log[-u]-Log[mu2],Log[-s/mu2]->Log[-s]-Log[mu2]};
RuleBox2={Log[-s]->Log[s/u]+Log[-u]};

Timing[
AppendTo[MasInt2,j[Base,1,0,1,0]->SolvejDGLm3[j[Base,1,0,1,0],s,nmaximal,2,GetBoundaries2[j[Base,1,0,1,0],s,2,Boundaries]]];
AppendTo[MasInt2,j[Base,0,1,0,1]->SolvejDGLm3[j[Base,0,1,0,1],u,nmaximal,2,GetBoundaries2[j[Base,0,1,0,1],u,2,Boundaries]]];
AppendTo[MasInt2,j[Base,1,1,1,0]->SolvejDGLm3[j[Base,1,1,1,0],s,nmaximal,2,GetBoundaries3[j[Base,1,1,1,0],s,2,Boundaries]]];
AppendTo[MasInt2,j[Base,1,0,1,1]->SolvejDGLm3[j[Base,1,0,1,1],s,nmaximal,2,GetBoundaries3[j[Base,1,0,1,1],s,2,Boundaries]]];
temp=SolvejDGLm3[j[Base,0,1,1,1],u,nmaximal,2,GetBoundaries3[j[Base,0,1,1,1],u,2,Boundaries]];
AppendTo[MasInt2,j[Base,0,1,1,1]->temp];
AppendTo[MasInt2,j[Base,1,1,0,1]->temp];
AppendTo[MasInt2,j[Base,1,1,1,1]->SolvejDGLm3[j[Base,1,1,1,1],u,nmaximal,0,GetBoundaries4[Boundaries]]];
MasInt2[[2]]=MasInt2[[2]]/.Rule2s;
MasInt2[[3]]=MasInt2[[3]]/.Rule2u;
MasInt2[[4]]=MasInt2[[4]]/.Rule3s;
MasInt2[[5]]=MasInt2[[5]]/.Rule3s;
MasInt2[[6]]=MasInt2[[6]]/.Rule3u;
MasInt2[[7]]=MasInt2[[7]]/.Rule2u;
MasInt2[[8]]=MasInt2[[8]]/.RuleBox/.RuleBox2;
Export["MasInt2FFF.m",MasInt2]
]












