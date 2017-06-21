(* ::Package:: *)

(*Developer:: adetokunbo adedoyin*)
(*contact:: adedoyin.adetokunbo@gmail.com*)
(*material_model:: bcj*)
(*model_dimension:: 3D *)

(*Parameter are for 304L Steel*)
(*BCJ-Model - Simple Shear*)
(*BCJ Model*)

\[Epsilon] = $MachineEpsilon;

(*material parameters*)
NU = 0.3; 
E1 = 206804;

(*yield parameters*)
f = 0.001; y = 100; v = 71.097195646; 
(*Isotropic Hardening*)
h = 200; rs = 0.00250746; rd = 0.001; 
(*Kinematic Hardening*)
H = 1300; Rs = 0.001; Rd = 0.001; 
(*material parameters definitions*)
MU = E1/(2*(1+NU)); 
LAMDA = E1*NU*(1 /( (1+NU)*(1-(2*NU)) ));

(*memory declarations*)
F[t_] := {{1, 1*t, 0},{0, 1, 0},{0, 0, 1}};
L[t_] := F'[t].Inverse[F[t]];
DT[t_] := (1/2)*(L[t] + Transpose[L[t]]);
WT[t_] := (1/2)*(L[t] - Transpose[L[t]]);
Sigma[t_] = Table[Subscript[sigma,i,j][t], {i,1,3},{j,1,3}];
Alpha[t_] = Table[Subscript[alpha,i,j][t], {i,1,3},{j,1,3}];
DP[t_] = Table[Subscript[dp,i,j][t], {i,1,3},{j,1,3}];

(*Tensoral definitions*)
NORMDP[t_] := Sqrt[2/3]*Sqrt[Tr[(DP[t]).Transpose[DP[t]]]];
NORMALPHA[t_] := Sqrt[2/3]*Sqrt[Tr[(Alpha[t]).Transpose[Alpha[t]]]]
Deviatoric[x_] := x-(1/3)*Tr[x]*IdentityMatrix[3];
NORMPHI[t_] := Sqrt[2/3]*Sqrt[3*Tr[(Sigma[t]-(2/3)*Alpha[t]).Transpose[(Sigma[t]-(2/3)*Alpha[t])]]];

(*DP*)
DP[t_]:= f*Sinh[(NORMPHI[t]-y-Kappa[t])/v]*((Sigma[t]-(2/3)*Alpha[t])/NORMPHI[t]);


(*Equations*)
eqns1 = Thread[ D[Flatten[Sigma[t]],t] == LAMDA*Tr[Deviatoric[DT[t]]-DP[t]]*Flatten[IdentityMatrix[3]] + 2*MU*Flatten[(Deviatoric[DT[t]]-DP[t])] + Flatten[WT[t].Sigma[t]] - Flatten[Sigma[t].WT[t]] ];
eqns2 = Thread[ D[Flatten[Alpha[t]],t] == h*Flatten[DP[t]] - ( rd*(2/3)*NORMALPHA[t] + rs )*(2/3)*Flatten[Alpha[t]]*(NORMALPHA[t]) + Flatten[WT[t].Alpha[t]] - Flatten[Alpha[t].WT[t]] ];
eqns3 = Thread[ D[Kappa[t],t] == h*NORMDP[t] - (Rd*(2/3)*NORMDP[t] + Rs)*(2/3)*(Kappa[t]^2) ];
eqns = Join[eqns1,eqns2,{eqns3}];

initc1=Thread[Flatten[Sigma[0]] == Flatten[{{0.0001,0.0001,0},{0.0001,0.0001,0},{0,0,0}}]];
initc2=Thread[Flatten[Alpha[0]] == Flatten[{{0.001,0,0},{0.001,0,0},{0,0,0}}]];
initc3 = Thread[ Kappa[0] == 0.009 ];
initc = Join[initc1,initc2,{initc3}];

func1 = Thread[ Flatten[ Sigma[t]] ];
func2 = Thread[ Flatten[ Alpha[t]] ];
func3 = Thread[ Kappa[t] ] ;
func = Join[func1,func2,{func3}];


(*Plotting Routine*)
intialT=0;
finalT=1;
(*Ppl*)
s = NDSolve[{eqns,initc},func,{t,intialT,finalT}];
List[{
Table[Plot[Evaluate[Subscript[sigma, 1,1][t]/.s],{t,intialT,finalT},PlotLabel->f[t],ImageSize->Scaled[0.2], PlotRange->All],{f,{Subscript[sigma, 1,1]}}],

Table[Plot[Evaluate[Subscript[sigma, 1,2][t]/.s],{t,intialT,finalT},PlotLabel->f[t],
ImageSize->Scaled[0.2],
AxesOrigin->{0,0},
GridLines->Automatic,
PlotStyle->{Black,Thick},PlotRange->All],
{f,{Subscript[sigma, 1,2]}}],

Table[Plot[Evaluate[Subscript[sigma, 1,3][t]/.s],{t,intialT,finalT},PlotLabel->f[t],ImageSize->Scaled[0.2], PlotRange->All],{f,{Subscript[sigma, 1,3]}}]}]

List[{
Table[Plot[Evaluate[Subscript[sigma, 2,1][t]/.s],{t,intialT,finalT},PlotLabel->f[t],
ImageSize->Scaled[0.2],
AxesOrigin->{0,0},
GridLines->Automatic,
PlotStyle->{Black,Thick},PlotRange->All],{f,{Subscript[sigma, 1,2]}}],

Table[Plot[Evaluate[Subscript[sigma, 2,2][t]/.s],{t,intialT,finalT},PlotLabel->f[t],ImageSize->Scaled[0.2], PlotRange->All],{f,{Subscript[sigma, 2,2]}}],
Table[Plot[Evaluate[Subscript[sigma, 2,3][t]/.s],{t,intialT,finalT},PlotLabel->f[t],ImageSize->Scaled[0.2], PlotRange->All],{f,{Subscript[sigma, 2,3]}}]
}]

List[{
Table[Plot[Evaluate[Subscript[sigma, 3,1][t]/.s],{t,intialT,finalT},PlotLabel->f[t],ImageSize->Scaled[0.2], PlotRange->All],{f,{Subscript[sigma, 3,1]}}],
Table[Plot[Evaluate[Subscript[sigma, 3,2][t]/.s],{t,intialT,finalT},PlotLabel->f[t],ImageSize->Scaled[0.2],PlotRange->All],{f,{Subscript[sigma, 3,2]}}],
Table[Plot[Evaluate[Subscript[sigma, 3,3][t]/.s],{t,intialT,finalT},PlotLabel->f[t],ImageSize->Scaled[0.2],PlotRange->All],{f,{Subscript[sigma, 3,3]}}]
}]




