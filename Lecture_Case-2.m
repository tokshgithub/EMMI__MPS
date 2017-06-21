(* ::Package:: *)

(*Developer:: adetokunbo adedoyin*)
(*contact:: adedoyin.adetokunbo@gmail.com*)
(*material_model:: bcj*)
(*model_dimension:: 1D *)

(*BCJ Model- With Isotropic Hardening (Case 2 )*)
\[Epsilon] = $MachineEpsilon; 

(*material parameters*)
E1 = 200000;
(*yield parameters*)
y = 100; f = 0.001; v = 71.097195646; 
(*Isotropic Hardening*)
H = 1300; Rs = 0.1; Rd = 0.001;


(*Auxillary Equations*)
(*plastic flow*)
DEpsillonP[t_]:= f*Sinh[(Abs[Sigma[t]-Alpha[t]]-Kappa[t]-y)/v]*Sign[Sigma[t]-Alpha[t]];


eqns1=Thread[\!\(
\*SubscriptBox[\(\[PartialD]\), \({t}\)]\(Epsillon[t]\)\)==StrainRate];
eqns2=Thread[\!\(
\*SubscriptBox[\(\[PartialD]\), \({t}\)]\(Sigma[t]\)\)==E1*(\!\(
\*SubscriptBox[\(\[PartialD]\), \({t}\)]\(Epsillon[t]\)\)-DEpsillonP[t])];
eqns3=Thread[\!\(
\*SubscriptBox[\(\[PartialD]\), \({t}\)]\(Alpha[t]\)\)==0.0];
eqns4=Thread[\!\(
\*SubscriptBox[\(\[PartialD]\), \({t}\)]\(Kappa[t]\)\)==H*DEpsillonP[t]-Rd Kappa[t]^2*DEpsillonP[t]-Rs Kappa[t]^2];
eqns=Join[{eqns1},{eqns2},{eqns3},{eqns4}];

initc1=Thread[Epsillon[\[Epsilon]]-0.00001==0];
initc2=Thread[Sigma[\[Epsilon]]==0.0001];
initc3=Thread[Alpha[\[Epsilon]]==0];
initc4=Thread[Kappa[\[Epsilon]]==0.0000001];
initc=Join[{initc1},{initc2},{initc3},{initc4}];

func1=Thread[Epsillon[t]];
func2=Thread[Sigma[t]];
func3=Thread[Alpha[t]];
func4=Thread[Kappa[t]];
func=Join[{func1},{func2},{func3},{func4}];

intialT=0;
finalT=1;
BCJ[StrainRate_?NumericQ] := NDSolve[{Join[{Thread[\!\(
\*SubscriptBox[\(\[PartialD]\), \(t\)]\(Epsillon[t]\)\) == StrainRate ]},{eqns2},{eqns3},{eqns4}],initc},func,{t,intialT,finalT},SolveDelayed->True][[1, 2, 2]]

P1=Plot[Evaluate[{BCJ[0.1],BCJ[0.3],BCJ[0.9]}],{t, intialT, finalT},
	PlotLabel->Style[Framed["Stress vs. Strain"],16,Black,Background->White],
	AxesLabel->{t,Sigma},
	ImageSize->Scaled[0.4],
	AxesOrigin->{0,0},
	GridLines->Automatic,
	PlotStyle->{{Red,Thick},{Blue,Thick},{Black,Thick}}];

intialT=0;finalT=0.3;
P2=Plot[Evaluate[{BCJ[1],BCJ[3],BCJ[9]}],{t, intialT, finalT},
	PlotLabel->Style[Framed["Stress vs. Strain"],16,Black,Background->White],
	AxesLabel->{t,Sigma},
	ImageSize->Scaled[0.4],
	AxesOrigin->{0,0},
	GridLines->Automatic,
	PlotStyle->{{Red,Thick},{Blue,Thick},{Black,Thick}}];
List[P1,"\t",P2]



