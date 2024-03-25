#!/usr/bin/env wolframscript

(* ::Section:: *)
(* Initialisation *)


(* ::Input:: *)
ClearAll;
Needs["ErrorBarPlots`"];
Needs["PlotLegends`"];


(* ::Input:: *)
NT := _SED_NT_ (*Temporal width*)
NS := _SED_NS_ (*Spatial width*)
beta := _SED_BETA_ (*Beta value*)
Digitsb := _SED_DIGITSB_ (*Number of digits in beta value*)
DPb := _SED_DPB_ (*Number of decimal points of beta value*)
Digitsm := _SED_DIGITSM_ (*Number of digits in bare mass*)
DPm := _SED_DPM_ (*Number of decimal points of bare mass*)
m := _SED_M_  (*Input mass*)
Nc := _SED_NC_ (*Number of colours*)
Rep := _SED_REP_ (*Representation: F=Fundamental, S=Symmetric, AS=Antisymmetric*)
Basepath := "_SED_BASEPATH_";
Outputpath := "_SED_OUTPUTPATH_";
FilenameSlug := StringJoin["S",ToString[NS],"T",ToString[NT],"B", ToString[NumberForm[beta,{Digitsb,DPb}]],"_m", ToString[Rep], ToString[NumberForm[m,{Digitsm,DPm}]]];

Nboot := 200 (*Number of bootstrap iterations*)
channels:={"g5","id","g0","g1","g2","g3","g0g1","g0g2","g0g3","g0g5","g5g1","g5g2","g5g3","g0g5g1","g0g5g2","g0g5g3","g5_g0g5_re" };
Operators:={"\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)","\[DoubleStruckOne]","\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(1\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(2\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(3\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(1\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(2\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(3\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(1\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(2\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(3\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(1\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(2\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(3\)]\)","\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(_\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)_re" };
plaquette=Import[StringJoin[Basepath,"/plaquette_", FilenameSlug, ".dat"],"Table"];
SeedRandom[_SED_SEEDRANDOM_, Method->"MersenneTwister"];
RAWDATA=Import[StringJoin[Basepath,"/correlators_", FilenameSlug, ".dat"],"Table"];
Cnfg=Length[RAWDATA]/Length[channels];
RandInt=Table[Table[Round[RandomReal[{0,1}]*(Cnfg-1)+1],{i,Nboot}],{j,Cnfg}];
plaqboot=Table[Mean[Table[plaquette[[RandInt[[i,j]]]],{i,Cnfg}]],{j,Nboot}];
plaqbootavg=Mean[plaqboot][[1]];
Casimir:=If[StringMatchQ[ToString[Rep],"F"],(1+Nc)/4,If[StringMatchQ[ToString[Rep],"AS"],Nc/2,If[StringMatchQ[ToString[Rep],"S"],Nc/2+1]]];Subscript[\[CapitalDelta], \[CapitalSigma]]:=-12.82; Subscript[\[CapitalDelta], \[Gamma]5\[Gamma]\[Mu]]:=-3.0;Subscript[\[CapitalDelta], \[Gamma]\[Mu]]:=-7.75
ZV=1+Casimir*(Subscript[\[CapitalDelta], \[CapitalSigma]]+Subscript[\[CapitalDelta], \[Gamma]\[Mu]])*Nc/(beta*plaqbootavg*8\[Pi]^2);
ZA=1+Casimir*(Subscript[\[CapitalDelta], \[CapitalSigma]]+Subscript[\[CapitalDelta], \[Gamma]5\[Gamma]\[Mu]])*Nc/(beta*plaqbootavg*8\[Pi]^2);
GradientFlow=_SED_GRADIENTFLOW_;
ChiSqCrit[dof_,percent_]:=Quiet[(\[Chi]/.NSolve[CDF[ChiSquareDistribution[dof],\[Chi]]-percent==0][[1]])/dof];
(*Function that returns the REDUCED \[Chi]^2 value for a given number of degrees of freedom (dof) that corresponds to a given certainty (percent)*)
BinWidth:=0.05
If[FileExistsQ[StringJoin[Outputpath,"/output_pseudoscalar.txt"]],DeleteFile[StringJoin[Outputpath,"/output_pseudoscalar.txt"]]];If[FileExistsQ[StringJoin[Outputpath,"/output_vector.txt"]],DeleteFile[StringJoin[Outputpath,"/output_vector.txt"]]];If[FileExistsQ[StringJoin[Outputpath,"/output_axialvector.txt"]],DeleteFile[StringJoin[Outputpath,"/output_axialvector.txt"]]];If[FileExistsQ[StringJoin[Outputpath,"/output_scalar.txt"]],DeleteFile[StringJoin[Outputpath,"/output_scalar.txt"]]];If[FileExistsQ[StringJoin[Outputpath,"/output_tensor.txt"]],DeleteFile[StringJoin[Outputpath,"/output_tensor.txt"]]];If[FileExistsQ[StringJoin[Outputpath,"/output_axialtensor.txt"]],DeleteFile[StringJoin[Outputpath,"/output_axialtensor.txt"]]];
If[FileExistsQ[StringJoin[Outputpath,"/pseudoscalarCSD.pdf"]],DeleteFile[StringJoin[Outputpath,"/pseudoscalarCSD.pdf"]]];If[FileExistsQ[StringJoin[Outputpath,"/vectorCSD.pdf"]],DeleteFile[StringJoin[Outputpath,"/vectorCSD.pdf"]]];If[FileExistsQ[StringJoin[Outputpath,"/axialvectorCSD.pdf"]],DeleteFile[StringJoin[Outputpath,"/axialvectorCSD.pdf"]]];If[FileExistsQ[StringJoin[Outputpath,"/scalarCSD.pdf"]],DeleteFile[StringJoin[Outputpath,"/scalarCSD.pdf"]]];If[FileExistsQ[StringJoin[Outputpath,"/tensorCSD.pdf"]],DeleteFile[StringJoin[Outputpath,"/tensorCSD.pdf"]]];If[FileExistsQ[StringJoin[Outputpath,"/axialtensorCSD.pdf"]],DeleteFile[StringJoin[Outputpath,"/axialtensorCSD.pdf"]]];
If[FileExistsQ[StringJoin[Outputpath, "/", FilenameSlug, "_pseudoscalar_masses_boots.csv"]],DeleteFile[StringJoin[Outputpath,"/",FilenameSlug,"_pseudoscalar_masses_boots.csv"]]];
If[FileExistsQ[StringJoin[Outputpath, "/", FilenameSlug, "_pseudoscalar_decayconsts_boots.csv"]],DeleteFile[StringJoin[Outputpath,"/",FilenameSlug,"_pseudoscalar_decayconsts_boots.csv"]]];
If[FileExistsQ[StringJoin[Outputpath, "/", FilenameSlug,"_vector_masses_boots.csv"]],DeleteFile[StringJoin[Outputpath, "/", FilenameSlug,"_vector_mass_boots.csv"]]];
If[FileExistsQ[StringJoin[Outputpath, "/", FilenameSlug,"_vector_decayconsts_boots.csv"]],DeleteFile[StringJoin[Outputpath, "/", FilenameSlug,"_vector_decayconsts_boots.csv"]]];
If[FileExistsQ[StringJoin[Outputpath, "/", FilenameSlug,"_axialvector_masses_boots.csv"]],DeleteFile[StringJoin[Outputpath, "/", FilenameSlug,"_axialvector_masses_boots.csv"]]];
If[FileExistsQ[StringJoin[Outputpath, "/", FilenameSlug,"_axialvector_decayconsts_boots.csv"]],DeleteFile[StringJoin[Outputpath, "/", FilenameSlug,"_axialvector_decayboots.csv"]]];
If[FileExistsQ[StringJoin[Outputpath, "/", FilenameSlug,"_scalar_masses_boots.csv"]],DeleteFile[StringJoin[Outputpath, "/", FilenameSlug,"_scalar_masses_boots.csv"]]];
If[FileExistsQ[StringJoin[Outputpath, "/", FilenameSlug,"_tensor_masses_boots.csv"]],DeleteFile[StringJoin[Outputpath, "/", FilenameSlug,"_tensor_masses_boots.csv"]]];
If[FileExistsQ[StringJoin[Outputpath, "/", FilenameSlug,"_axialtensor_masses_boots.csv"]],DeleteFile[StringJoin[Outputpath, "/", FilenameSlug,"_axialtensor_masses_boots.csv"]]];
(*Remove old data files to prevent them being used in subsequent plots*)


(* ::Section:: *)
(* Pseudoscalar *)


(* ::Input:: *)
UsePS=_SED_USEPS_;
AutomaticPS=_SED_AUTOMATICPS_;
MinimumPlateau=3;


(* ::Input:: *)
ch=1;
g5import=Table[RAWDATA[[Length[channels]*(i-1)+ch]],{i,1,Cnfg}];
For[i=1,i<=Cnfg,i++,g5import[[i]]=Append[g5import[[i]],g5import[[i,3]]]];


(* ::Input:: *)
ch=17;
g5g0g5import=Table[RAWDATA[[Length[channels]*(i-1)+ch]],{i,1,Cnfg}];
For[i=1,i<=Cnfg,i++,g5g0g5import[[i]]=Append[g5g0g5import[[i]],g5g0g5import[[i,3]]]];


(* ::Input:: *)
corrg5=Table[Table[NS^3(g5import[[j,i+2]]+g5import[[j,NT+4-i]])/2,{i,1,NT/2+1}],{j,1,Cnfg}];
corrg5g0g5=Table[Table[NS^3(g5g0g5import[[j,i+2]]-g5g0g5import[[j,NT+4-i]])/2,{i,1,NT/2+1}],{j,1,Cnfg}];
corrg5g0g5[[1,1]]=corrg5g0g5[[NT/2+1,NT/2+1]]=1;(*The top-left and bottom-right terms are set to one to allow for the computation of the inverse*)


(* ::Input:: *)
corrg5=Transpose[corrg5];
corrg5g0g5=Transpose[corrg5g0g5];


(* ::Input:: *)
RandInt=Table[Table[Round[RandomReal[{0,1}]*(Cnfg-1)+1],{i,Nboot}],{j,Cnfg}];
corrg5boot=Table[Table[Mean[Table[corrg5[[k,RandInt[[i,j]]]],{i,Cnfg}]],{j,Nboot}],{k,NT/2+1}];
corrg5bootavg=Table[{k-1,Mean[corrg5boot[[k]]],StandardDeviation[corrg5boot[[k]]]},{k,NT/2+1}];


(* ::Input:: *)
RandInt=Table[Table[Round[RandomReal[{0,1}]*(Cnfg-1)+1],{i,Nboot}],{j,Cnfg}];
corrg5g0g5boot=Table[Table[Mean[Table[corrg5g0g5[[k,RandInt[[i,j]]]],{i,Cnfg}]],{j,Nboot}],{k,NT/2+1}];
corrg5g0g5bootavg=Table[{k-1,Mean[corrg5g0g5boot[[k]]],StandardDeviation[corrg5g0g5boot[[k]]]},{k,NT/2+1}];


(* ::Input:: *)
meffg5=Table[{i,Mean[Table[ArcCosh[(corrg5boot[[i,j]]+corrg5boot[[i+2,j]])/(2*corrg5boot[[i+1,j]])],{j,Nboot}]],StandardDeviation[Table[ArcCosh[(corrg5boot[[i,j]]+corrg5boot[[i+2,j]])/(2*corrg5boot[[i+1,j]])],{j,Nboot}]]},{i,1,NT/2-1}];


(* ::Input:: *)
g5meffplot=ErrorListPlot[Table[{{meffg5[[t,1]],meffg5[[t,2]]},ErrorBar[meffg5[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->Automatic, FrameLabel->{"t","\!\(\*SubsuperscriptBox[\(am\), \(eff\), \(PS\)]\)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->Automatic,AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[0,0,1],PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]]]]


(* ::Input:: *)
Export[StringJoin[Outputpath,"/pseudoscalar.pdf"],g5meffplot,"PDF"];


(* ::Input:: *)
g5corrplot=ErrorListPlot[Table[{{corrg5bootavg[[t,1]],corrg5bootavg[[t,2]]},ErrorBar[corrg5bootavg[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->{"o",15}, FrameLabel->{"t","\!\(\*SubscriptBox[\(C\), \(\[Gamma]5\)]\)(t)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->{{0,NT},All},AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[1,0,0]]


(* ::Input:: *)
g5g0g5corrplot=ErrorListPlot[Table[{{corrg5g0g5bootavg[[t,1]],corrg5g0g5bootavg[[t,2]]},ErrorBar[corrg5g0g5bootavg[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->{"o",15}, FrameLabel->{"t","\!\(\*SubscriptBox[\(C\), \(\[Gamma]5_\[Gamma]0\[Gamma]5\)]\)(t)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->{{0,NT},All},AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[1,0,0]]


(* ::Input:: *)
Clear[f,g]
f[a_,b_,c_,t_]:=(a^2/b)*(Exp[-b*t]+Exp[-b(NT-t)]);
g[a_,b_,c_,t_]:=(a*c)*(Exp[-b*t]-Exp[-b(NT-t)]);
(*Fitting functions for matrix element (a), mass (b) and decay constant (c)*)
BootstrapPS[IPps_,FPps_]:=Table[w1=i;tmp = NMinimize[Sum[Sum[(corrg5boot[[j,i]]-f[a,b,c,corrg5bootavg[[j,1]]])*invg5cov[[j,k]]*(corrg5boot[[k,i]]-f[a,b,c,corrg5bootavg[[k,1]]])+(corrg5g0g5boot[[j,i]]-g[a,b,c,corrg5g0g5bootavg[[j,1]]])*invg5g0g5cov[[j,k]]*(corrg5g0g5boot[[k,i]]-g[a,b,c,corrg5g0g5bootavg[[k,1]]]), {j, IPps, FPps}], {k, IPps, FPps}], {a,b,c}, Method -> "NelderMead"];
dof=2(FPps-IPps+1)-3;
 {a/. tmp[[2,1]], b /. tmp[[2,2]],c/.tmp[[2,3]],tmp[[1]]/dof}, {i, 1, Nboot}];
(*Outputs {{a,b,c,\[Chi]^2/dof},...} Nboot times*)
BootErrorsPS[input_]:=Table[DistributionFitTest[Abs[input[[;;,i]]],Automatic,"HypothesisTestData"]["FittedDistribution"][[2]],{i,1,Dimensions[input][[2]]}]
CurveFitPS[IPps_,FPps_]:={tmp = NMinimize[Sum[Sum[(corrg5bootavg[[j,2]]-f[a,b,c,corrg5bootavg[[j,1]]])*invg5cov[[j,k]]*(corrg5bootavg[[k,2]]-f[a,b,c,corrg5bootavg[[k,1]]])+(corrg5g0g5bootavg[[j,2]]-g[a,b,c,corrg5g0g5bootavg[[j,1]]])*invg5g0g5cov[[j,k]]*(corrg5g0g5bootavg[[k,2]]-g[a,b,c,corrg5g0g5bootavg[[k,1]]]), {j, IPps, FPps}], {k, IPps, FPps}], {a,b,c}, Method -> "NelderMead"];
dof=2(FPps-IPps+1)-3;
Abs[a]/. tmp[[2,1]], Abs[b] /. tmp[[2,2]],Abs[c]/.tmp[[2,3]],tmp[[1]]/dof};
(*Outputs {|a|,|b|,|c|,\[Chi]^2/dof}*)


(* ::Input:: *)
g5covariance=Table[Table[Mean[Table[(corrg5boot[[j,i]]-corrg5bootavg[[j,2]])(corrg5boot[[k,i]]-corrg5bootavg[[k,2]]),{i,Nboot}]],{j,NT/2+1}],{k,NT/2+1}];
g5g0g5covariance=Table[Table[Mean[Table[(corrg5g0g5boot[[j,i]]-corrg5g0g5bootavg[[j,2]])(corrg5g0g5boot[[k,i]]-corrg5g0g5bootavg[[k,2]]),{i,Nboot}]],{j,NT/2+1}],{k,NT/2+1}];
invg5cov=Inverse[g5covariance];
invg5g0g5cov=Inverse[g5g0g5covariance];


(* ::Input:: *)
(*Find the longest sequence of t values for which the effective mass is both real and positive.*)
tmpMasses={};
For[i=1,i<=Dimensions[meffg5][[1]],++i,If[Element[meffg5[[i,2]],Reals],AppendTo[tmpMasses,meffg5[[i,;;]]]]]
setsplit=Split[tmpMasses[[;;,1]],#2-#1==1&];
lengths=Table[Length[setsplit[[i,;;]]],{i,Dimensions[setsplit][[1]]}];
place=Position[lengths,Max[lengths]][[1,1]];
window=setsplit[[place]]
(*If there are fewer than a specified number of points that satisfy the above, the fitting must be done by hand.*)
If[Length[window]<=0.4*Length[meffg5],{Print["Plateau fitting must be done manually!"],AutomaticPS=False;}]


(* ::Input:: *)
If[AutomaticPS&&UsePS,TInitial=window[[1]];TFinal=window[[-1]];
tmp=CurveFitPS[TInitial,TFinal];
chisq=tmp[[-1]];
For[t=TInitial+1,t<TFinal,++t,tmp=CurveFitPS[t,TFinal];dofps=2(TFinal-t+1)-3;
If[tmp[[-1]]<chisq&&tmp[[-1]]>ChiSqCrit[dofps,0.90],chisq=tmp[[-1]],{TInitial=t-1,Break[]}]];TInitial=t-1;]


(* ::Input:: *)
If[AutomaticPS&&UsePS,{For[cert=0.90,cert<=0.95,cert+=0.05,{If[AutomaticPS&&UsePS,Variation=2;Tests={};
For[y=TFinal,y>=TFinal-Variation,--y,For[x=Max[1,TInitial-Variation],x<=y-MinimumPlateau,++x,
tmp=CurveFitPS[x,y];
dofps=2(y-x+1)-3;
chisq=tmp[[-1]];
If[chisq>ChiSqCrit[dofps,cert],Continue[]];
chisq/=(y-x+1)^2;
AppendTo[Tests,{chisq,x,y}]]];
If[Length[Tests]!=0,{FitRange=Tests[[Position[Tests[[;;,1]],Min[Tests[[;;,1]]]][[1,1]]]],Break[]},Print["No valid fit ranges at cert = "<>ToString[cert]]];]}]}];
(*Select the fit range for which \[Chi]^2/dof(tf-ti+1)^2 is minimised. This attempts to simultaneously maximise the legnth of the fit range and minimise \[Chi]^2/dof. In some cases the initial certainty may be too restrictive. If so, the certainty is changed to widen the search.*)


(* ::Input:: *)
If[Length[Tests]==0,AutomaticPS=False];(*If cert is still too restrictive, the results have to be done by hand.*)


(* ::Input:: *)
If[AutomaticPS&&UsePS,{{IPps=FitRange[[2]],FPps=FitRange[[3]]},dofps=2(FPps-IPps+1)-3}]
If[Not[AutomaticPS]&&UsePS,{{IPps=_SED_IPPS_,FPps=_SED_FPPS_},dofps=2(FPps-IPps+1)-3}]


(* ::Input:: *)
ProgressIndicator[Dynamic[w1],{1,Nboot}]


(* ::Input:: *)
If[UsePS,Provisional=CurveFitPS[IPps,FPps]](*These are the provisional results that are subject to change during the robustness check.*)
(*Matrix element, mass, decay-constant, \[Chi]^2/dof*)


(* ::Input:: *)
If[UsePS,BootsProvisional=BootstrapPS[IPps,FPps]];(*BootsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[UsePS,ErrsProvisional=BootErrorsPS[BootsProvisional]](*ErrsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[AutomaticPS&&UsePS,For[t=IPps+1,t<FPps+2-MinimumPlateau,++t,{temp=CurveFitPS[t,FPps];tempBoots=BootstrapPS[t,FPps];tempErr=BootErrorsPS[tempBoots];If[Abs[temp[[2]]-Provisional[[2]]]<ErrsProvisional[[2]],{PSFinal=Provisional,PSBoots=BootsProvisional,PSErr=ErrsProvisional,Break[]},{Provisional=temp,BootsProvisional=tempBoots,ErrsProvisional=tempErr}]}]]


(* ::Input:: *)
If[Not[AutomaticPS]&&UsePS,{PSFinal=Provisional,PSBoots=BootsProvisional,PSErr=ErrsProvisional}];


(* ::Input:: *)
{{IPps,FPps},dofps}


(* ::Input:: *)
If[UsePS,FINALSPS=ArrayReshape[Riffle[PSFinal,PSErr],{Length[PSFinal],2}]]
If[UsePS,FinalChiSq=FINALSPS[[-1,1]]];


(* ::Input:: *)
If[UsePS,For[i=1,i<=Nboot,++i,Write[StringJoin[Outputpath, "/", FilenameSlug,"_pseudoscalar_masses_boots.csv"],(PSBoots[[i,2]]*GradientFlow)^2]], Write[StringJoin[Outputpath, "/", FilenameSlug,"_pseudoscalar_masses_boots.csv"], ""]];
If[UsePS,For[i=1,i<=Nboot,++i,Write[StringJoin[Outputpath, "/", FilenameSlug,"_pseudoscalar_decayconsts_boots.csv"],(ZA*PSBoots[[i,3]]*GradientFlow)^2]], Write[StringJoin[Outputpath, "/", FilenameSlug,"_pseudoscalar_decayconsts_boots.csv"], ""]];
(*Write all bootstraps to an output file for use in the chiral extrapolation.*)


(* ::Input:: *)
If[UsePS,i=1;data=Abs[PSBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[PSFinal[[i]],PSErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[PSFinal[[i]]],", ",ToString[PSErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*TemplateBox[{\"0\"},\n\"Bra\"]\)\!\(\*SubscriptBox[\(\[ScriptCapitalO]\), \(PS\)]\)\!\(\*TemplateBox[{\"PS\"},\n\"Ket\"]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]
If[UsePS,i=2;data=Abs[PSBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[PSFinal[[i]],PSErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[PSFinal[[i]]],", ",ToString[PSErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*SubscriptBox[\(m\), \(PS\)]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]
If[UsePS,i=3;data=Abs[PSBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[PSFinal[[i]],PSErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[PSFinal[[i]]],", ",ToString[PSErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*SubscriptBox[\(f\), \(PS\)]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]


(* ::Input:: *)
If[UsePS,FINALPLOTPS=Plot[{FINALSPS[[2,1]]+FINALSPS[[2,2]],FINALSPS[[2,1]]-FINALSPS[[2,2]]},{t,IPps-0.5,FPps+0.5},Filling->True];
If[AutomaticPS,FullPlotPS=Show[g5meffplot,FINALPLOTPS,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]]]]]]
If[Not[AutomaticPS],FullPlotPS=Show[g5meffplot,FINALPLOTPS,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]],"\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"Manual\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\")\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)"]]]


(* ::Input:: *)
If[UsePS,Export[StringJoin[Outputpath,"/pseudoscalar.pdf"],FullPlotPS,"PDF"],Export[StringJoin[Outputpath,"/pseudoscalar.pdf"], ""]];
If[UsePS,HistogramPlotPS=Histogram[PSBoots[[;;,-1]],{BinWidth},"Probability",ImageSize->600,ChartStyle->Darker[Blue],FrameLabel->{"\!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.","Normalised Frequency"},LabelStyle->Directive[20,Black],Frame->True,PlotRange->Automatic,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT], ", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=",ToString[NumberForm[m,{Digitsm,DPm}]]]]]


(* ::Input:: *)
If[UsePS,Export[StringJoin[Outputpath,"/pseudoscalarCSD.pdf"],HistogramPlotPS,"PDF"], Export[StringJoin[Outputpath,"/pseudoscalarCSD.pdf"], ""]];


(* ::Input:: *)
(*Square measurements, if necessary, accounting for propagation of errors and renormalise decay constants.*)
If[UsePS,amplitudePS=Around[FINALSPS[[1,1]],FINALSPS[[1,2]]];
massPS=Around[FINALSPS[[2,1]],FINALSPS[[2,2]]];
fPS=Around[FINALSPS[[3,1]],FINALSPS[[3,2]]];
amplitudePS=amplitudePS^2;
massPS=(GradientFlow*massPS)^2;
fPS=(GradientFlow*fPS)^2;]


(* ::Input:: *)
If[UsePS,FINALSPS[[1,1]]=amplitudePS["Value"];
FINALSPS[[1,2]]=amplitudePS["Error"];
FINALSPS[[2,1]]=massPS["Value"];
FINALSPS[[2,2]]=massPS["Error"];
FINALSPS[[3,1]]=ZA^2*fPS["Value"];
FINALSPS[[3,2]]=ZA^2*fPS["Error"];]


(* ::Input:: *)
If[UsePS,Export[StringJoin[Outputpath,"/output_pseudoscalar.txt"],{Flatten[FINALSPS],ToString[IPps],ToString[FPps]},"Table","FieldSeparators"->"\n"], Export[StringJoin[Outputpath,"/output_pseudoscalar.txt"], ""]];


(* ::Section:: *)
(* Vector *)


(* ::Input:: *)
UseV=_SED_USEV_;
AutomaticV=_SED_AUTOMATICV_;
MinimumPlateau=3;


(* ::Input:: *)
ch=4;
g1import=Table[RAWDATA[[Length[channels]*(i-1)+ch]],{i,1,Cnfg}];
g2import=Table[RAWDATA[[Length[channels]*(i-1)+ch+1]],{i,1,Cnfg}];
g3import=Table[RAWDATA[[Length[channels]*(i-1)+ch+2]],{i,1,Cnfg}];


(* ::Input:: *)
For[i=1,i<=Cnfg,i++,g1import[[i]]=Append[g1import[[i]],g1import[[i,3]]]];
For[i=1,i<=Cnfg,i++,g2import[[i]]=Append[g2import[[i]],g2import[[i,3]]]];
For[i=1,i<=Cnfg,i++,g3import[[i]]=Append[g3import[[i]],g3import[[i,3]]]];


(* ::Input:: *)
corrV=Table[Table[NS^3(g1import[[j,i+2]]+g2import[[j,i+2]]+g3import[[j,i+2]]+g1import[[j,NT+4-i]]+g2import[[j,NT+4-i]]+g3import[[j,NT+4-i]])/6,{i,1,NT/2+1}],{j,1,Cnfg}];


(* ::Input:: *)
corrV=Transpose[corrV];


(* ::Input:: *)
RandInt=Table[Table[Round[RandomReal[{0,1}]*(Cnfg-1)+1],{i,Nboot}],{j,Cnfg}];
corrVboot=Table[Table[Mean[Table[corrV[[k,RandInt[[i,j]]]],{i,Cnfg}]],{j,Nboot}],{k,NT/2+1}];
corrVbootavg=Table[{k-1,Mean[corrVboot[[k]]],StandardDeviation[corrVboot[[k]]]},{k,NT/2+1}];


(* ::Input:: *)
meffV=Table[{i,Mean[Table[ArcCosh[(corrVboot[[i,j]]+corrVboot[[i+2,j]])/(2*corrVboot[[i+1,j]])],{j,Nboot}]],StandardDeviation[Table[ArcCosh[(corrVboot[[i,j]]+corrVboot[[i+2,j]])/(2*corrVboot[[i+1,j]])],{j,Nboot}]]},{i,1,NT/2-1}];


(* ::Input:: *)
meffVplot=ErrorListPlot[Table[{{meffV[[t,1]],meffV[[t,2]]},ErrorBar[meffV[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->Automatic, FrameLabel->{"t","\!\(\*SubsuperscriptBox[\(am\), \(eff\), \(V\)]\)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->Automatic,AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[0,0,1],PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]]]]


(* ::Input:: *)
Export[StringJoin[Outputpath,"/vector.pdf"],meffVplot,"PDF"];


(* ::Input:: *)
Vcorrplot=ErrorListPlot[Table[{{corrVbootavg[[t,1]],corrVbootavg[[t,2]]},ErrorBar[corrVbootavg[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->{"o",15}, FrameLabel->{"t","\!\(\*SubscriptBox[\(C\), \(\[Gamma]k\)]\)(t)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->{{0,NT},All},AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[1,0,0]]


(* ::Input:: *)
Clear[f]
f[a_,b_,t_]:=(a*b)*(Exp[-b*t]+Exp[-b(NT-t)]);
(*Fitting function for the decay-constant squared (a) and the mass (b)*)
BootstrapV[IPv_,FPv_]:=Table[w2=i;tmp = NMinimize[Sum[Sum[(corrVboot[[j,i]]-f[a,b,corrVbootavg[[j,1]]])*invVcov[[j,k]]*(corrVboot[[k,i]]-f[a,b,corrVbootavg[[k,1]]]), {j, IPv, FPv}], {k, IPv, FPv}], {a,b}, Method -> "NelderMead"];
dof=(FPv-IPv+1)-2;
 {a/. tmp[[2,1]], b /. tmp[[2,2]],tmp[[1]]/dof}, {i, 1, Nboot}];
(*Outputs {{a,b,\[Chi]^2/dof},...} Nboot times*)
BootErrorsV[input_]:=Table[DistributionFitTest[Abs[input[[;;,i]]],Automatic,"HypothesisTestData"]["FittedDistribution"][[2]],{i,1,Dimensions[input][[2]]}]
CurveFitV[IPv_,FPv_]:={tmp = NMinimize[Sum[Sum[(corrVbootavg[[j,2]]-f[a,b,corrVbootavg[[j,1]]])*invVcov[[j,k]]*(corrVbootavg[[k,2]]-f[a,b,corrVbootavg[[k,1]]]), {j, IPv, FPv}], {k, IPv, FPv}], {a,b}, Method -> "NelderMead"];
dof=(FPv-IPv+1)-2;
Abs[a]/. tmp[[2,1]], Abs[b] /. tmp[[2,2]],tmp[[1]]/dof};
(*Outputs {|a|,|b|,\[Chi]^2/dof}*)


(* ::Input:: *)
Vcovariance=Table[Table[Mean[Table[(corrVboot[[j,i]]-corrVbootavg[[j,2]])(corrVboot[[k,i]]-corrVbootavg[[k,2]]),{i,Nboot}]],{j,NT/2+1}],{k,NT/2+1}];
invVcov=Inverse[Vcovariance];


(* ::Input:: *)
(*Find the longest sequence of t values for which the effective mass is both real and positive.*)
tmpMasses={};
For[i=1,i<=Dimensions[meffV][[1]],++i,If[Element[meffV[[i,2]],Reals],AppendTo[tmpMasses,meffV[[i,;;]]]]]
setsplit=Split[tmpMasses[[;;,1]],#2-#1==1&];
lengths=Table[Length[setsplit[[i,;;]]],{i,Dimensions[setsplit][[1]]}];
place=Position[lengths,Max[lengths]][[1,1]];
window=setsplit[[place]]
(*If there are fewer than a specified number of points that satisfy the above, the fitting must be done by hand.*)
If[Length[window]<=0.4*Length[meffV],{Print["Plateau fitting must be done manually!"],AutomaticV=False;}]


(* ::Input:: *)
If[AutomaticV&&UseV,TInitial=window[[1]];TFinal=window[[-1]];
tmp=CurveFitV[TInitial,TFinal];
chisq=tmp[[-1]];
For[t=TInitial+1,t<TFinal,++t,tmp=CurveFitV[t,TFinal];dofv=(TFinal-t+1)-2;
If[tmp[[-1]]<chisq&&tmp[[-1]]>ChiSqCrit[dofv,0.90],chisq=tmp[[-1]],{TInitial=t-1,Break[]}]];TInitial=t-1;]


(* ::Input:: *)
If[AutomaticV&&UseV,{For[cert=0.90,cert<=0.95,cert+=0.05,{If[AutomaticV&&UseV,Variation=2;Tests={};
For[y=TFinal,y>=TFinal-Variation,--y,For[x=Max[1,TInitial-Variation],x<=y-MinimumPlateau,++x,
tmp=CurveFitV[x,y];
dofv=(y-x+1)-2;
chisq=tmp[[-1]];
If[chisq>ChiSqCrit[dofv,cert],Continue[]];
chisq/=(y-x+1)^2;
AppendTo[Tests,{chisq,x,y}]]];
If[Length[Tests]!=0,{FitRange=Tests[[Position[Tests[[;;,1]],Min[Tests[[;;,1]]]][[1,1]]]],Break[]},Print["No valid fit ranges at cert = "<>ToString[cert]]];]}]}];
(*Select the fit range for which \[Chi]^2/dof(tf-ti+1)^2 is minimised. This attempts to simultaneously maximise the legnth of the fit range and minimise \[Chi]^2/dof. In some cases the initial certainty may be too restrictive. If so, the certainty is changed to widen the search.*)


(* ::Input:: *)
If[Length[Tests]==0,AutomaticV=False];(*If cert is still too restrictive, the results have to be done by hand.*)


(* ::Input:: *)
If[AutomaticV&&UseV,{{IPv=FitRange[[2]],FPv=FitRange[[3]]},dofv=(FPv-IPv+1)-2}]
If[Not[AutomaticV]&&UseV,{{IPv=_SED_IPV_,FPv=_SED_FPV_},dofv=(FPv-IPv+1)-2}]


(* ::Input:: *)
ProgressIndicator[Dynamic[w2],{1,Nboot}]


(* ::Input:: *)
If[UseV,Provisional=CurveFitV[IPv,FPv]](*These are the provisional results that are subject to change during the robustness check.*)
(*Squared decay-constant, mass, \[Chi]^2/dof*)


(* ::Input:: *)
If[UseV,BootsProvisional=BootstrapV[IPv,FPv]];(*BootsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[UseV,ErrsProvisional=BootErrorsV[BootsProvisional]](*ErrsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[AutomaticV&&UseV,For[t=IPv+1,t<FPv+2-MinimumPlateau,++t,{temp=CurveFitV[t,FPv];tempBoots=BootstrapV[t,FPv];tempErr=BootErrorsV[tempBoots];If[Abs[temp[[2]]-Provisional[[2]]]<ErrsProvisional[[2]],{VFinal=Provisional,VBoots=BootsProvisional,VErr=ErrsProvisional,Break[]},{Provisional=temp,BootsProvisional=tempBoots,ErrsProvisional=tempErr}]}]]


(* ::Input:: *)
If[Not[AutomaticV]&&UseV,{VFinal=Provisional,VBoots=BootsProvisional,VErr=ErrsProvisional}];


(* ::Input:: *)
{{IPv,FPv},dofv}


(* ::Input:: *)
If[UseV,FINALSV=ArrayReshape[Riffle[VFinal,VErr],{Length[VFinal],2}]]
If[UseV,FinalChiSq=FINALSV[[-1,1]]];


(* ::Input:: *)
If[UseV,For[i=1,i<=Nboot,++i,Write[StringJoin[Outputpath, "/", FilenameSlug,"_vector_masses_boots.csv"],(VBoots[[i,2]]*GradientFlow)^2]], Write[StringJoin[Outputpath, "/", FilenameSlug,"_vector_masses_boots.csv"], ""]];
If[UseV,For[i=1,i<=Nboot,++i,Write[StringJoin[Outputpath, "/", FilenameSlug,"_vector_decayconsts_boots.csv"],VBoots[[i,1]]*(ZV*GradientFlow)^2]], Write[StringJoin[Outputpath, "/", FilenameSlug,"_vector_decayconsts_boots.csv"], ""]];
(*Write all bootstraps to an output file for use in the chiral extrapolation.*)


(* ::Input:: *)
If[UseV,i=1;data=Abs[VBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[VFinal[[i]],VErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[VFinal[[i]]],", ",ToString[VErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*SubsuperscriptBox[\(f\), \(V\), \(2\)]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]
If[UseV,i=2;data=Abs[VBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[VFinal[[i]],VErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[VFinal[[i]]],", ",ToString[VErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*SubscriptBox[\(m\), \(V\)]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]


(* ::Input:: *)
If[UseV,FINALPLOTV=Plot[{FINALSV[[2,1]]+FINALSV[[2,2]],FINALSV[[2,1]]-FINALSV[[2,2]]},{t,IPv-0.5,FPv+0.5},Filling->True];
If[AutomaticV,FullPlotV=Show[meffVplot,FINALPLOTV,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]]]]]]
If[Not[AutomaticV],FullPlotV=Show[meffVplot,FINALPLOTV,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]],"\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"Manual\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\")\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)"]]]


(* ::Input:: *)
If[UseV,Export[StringJoin[Outputpath,"/vector.pdf"],FullPlotV,"PDF"], Export[StringJoin[Outputpath,"/vector.pdf"], ""]];
If[UseV,HistogramPlotV=Histogram[VBoots[[;;,-1]],{BinWidth},"Probability",ImageSize->600,ChartStyle->Darker[Blue],FrameLabel->{"\!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.","Normalised Frequency"},LabelStyle->Directive[20,Black],Frame->True,PlotRange->Automatic,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT], ", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=",ToString[NumberForm[m,{Digitsm,DPm}]]]]]


(* ::Input:: *)
If[UseV,Export[StringJoin[Outputpath,"/vectorCSD.pdf"],HistogramPlotV,"PDF"], Export[StringJoin[Outputpath,"/vectorCSD.pdf"], ""]];


(* ::Input:: *)
(*Square measurements, if necessary, accounting for propagation of errors and renormalise decay constants.*)
If[UseV,fV2=Around[FINALSV[[1,1]],FINALSV[[1,2]]];
massV=Around[FINALSV[[2,1]],FINALSV[[2,2]]];
fV2=fV2*(GradientFlow)^2;
massV=(GradientFlow*massV)^2;
FINALSV[[1,1]]=ZV^2*fV2["Value"];
FINALSV[[1,2]]=ZV^2*fV2["Error"];
FINALSV[[2,1]]=massV["Value"];
FINALSV[[2,2]]=massV["Error"];]


(* ::Input:: *)
If[UseV,Export[StringJoin[Outputpath,"/output_vector.txt"],{Flatten[FINALSV],ToString[IPv],ToString[FPv]},"Table","FieldSeparators"->"\n"], Export[StringJoin[Outputpath,"/output_vector.txt"], ""]];


(* ::Section:: *)
(* Axial-Vector *)


(* ::Input:: *)
UseAV=_SED_USEAV_;
AutomaticAV=_SED_AUTOMATICAV_;
MinimumPlateau=3;


(* ::Input:: *)
ch=11;
g5g1import=Table[RAWDATA[[Length[channels]*(i-1)+ch]],{i,1,Cnfg}];
g5g2import=Table[RAWDATA[[Length[channels]*(i-1)+ch+1]],{i,1,Cnfg}];
g5g3import=Table[RAWDATA[[Length[channels]*(i-1)+ch+2]],{i,1,Cnfg}];


(* ::Input:: *)
For[i=1,i<=Cnfg,i++,g5g1import[[i]]=Append[g5g1import[[i]],g5g1import[[i,3]]]];
For[i=1,i<=Cnfg,i++,g5g2import[[i]]=Append[g5g2import[[i]],g5g2import[[i,3]]]];
For[i=1,i<=Cnfg,i++,g5g3import[[i]]=Append[g5g3import[[i]],g5g3import[[i,3]]]];


(* ::Input:: *)
corrAV=Table[Table[NS^3(g5g1import[[j,i+2]]+g5g2import[[j,i+2]]+g5g3import[[j,i+2]]+g5g1import[[j,NT+4-i]]+g5g2import[[j,NT+4-i]]+g5g3import[[j,NT+4-i]])/6,{i,1,NT/2+1}],{j,1,Cnfg}];


(* ::Input:: *)
corrAV=Transpose[corrAV];


(* ::Input:: *)
RandInt=Table[Table[Round[RandomReal[{0,1}]*(Cnfg-1)+1],{i,Nboot}],{j,Cnfg}];
corrAVboot=Table[Table[Mean[Table[corrAV[[k,RandInt[[i,j]]]],{i,Cnfg}]],{j,Nboot}],{k,NT/2+1}];
corrAVbootavg=Table[{k-1,Mean[corrAVboot[[k]]],StandardDeviation[corrAVboot[[k]]]},{k,NT/2+1}];


(* ::Input:: *)
meffAV=Table[{i,Mean[Table[ArcCosh[(corrAVboot[[i,j]]+corrAVboot[[i+2,j]])/(2*corrAVboot[[i+1,j]])],{j,Nboot}]],StandardDeviation[Table[ArcCosh[(corrAVboot[[i,j]]+corrAVboot[[i+2,j]])/(2*corrAVboot[[i+1,j]])],{j,Nboot}]]},{i,1,NT/2-1}];


(* ::Input:: *)
meffAVplot=ErrorListPlot[Table[{{meffAV[[t,1]],meffAV[[t,2]]},ErrorBar[meffAV[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->Automatic, FrameLabel->{"t","\!\(\*SubsuperscriptBox[\(am\), \(eff\), \(AV\)]\)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->Automatic,AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[0,0,1],PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]]]]


(* ::Input:: *)
Export[StringJoin[Outputpath,"/axialvector.pdf"],meffAVplot,"PDF"];


(* ::Input:: *)
AVcorrplot=ErrorListPlot[Table[{{corrAVbootavg[[t,1]],corrAVbootavg[[t,2]]},ErrorBar[corrAVbootavg[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->{"o",15}, FrameLabel->{"t","\!\(\*SubscriptBox[\(C\), \(\[Gamma]5\[Gamma]k\)]\)(t)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->{{0,NT},All},AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[1,0,0]]


(* ::Input:: *)
Clear[f]
f[a_,b_,t_]:=(a*b)*(Exp[-b*t]+Exp[-b(NT-t)]);
(*Fitting function for the decay-constant squared (a) and the mass (b)*)
BootstrapAV[IPav_,FPav_]:=Table[w3=i;tmp = NMinimize[Sum[Sum[(corrAVboot[[j,i]]-f[a,b,corrAVbootavg[[j,1]]])*invAVcov[[j,k]]*(corrAVboot[[k,i]]-f[a,b,corrAVbootavg[[k,1]]]), {j, IPav, FPav}], {k, IPav, FPav}], {a,b}, Method -> "NelderMead"];
dof=(FPav-IPav+1)-2;
 {a/. tmp[[2,1]], b /. tmp[[2,2]],tmp[[1]]/dof}, {i, 1, Nboot}];
(*Outputs {{a,b,\[Chi]^2/dof},...} Nboot times*)
BootErrorsAV[input_]:=Table[DistributionFitTest[Abs[input[[;;,i]]],Automatic,"HypothesisTestData"]["FittedDistribution"][[2]],{i,1,Dimensions[input][[2]]}]
CurveFitAV[IPav_,FPav_]:={tmp = NMinimize[Sum[Sum[(corrAVbootavg[[j,2]]-f[a,b,corrAVbootavg[[j,1]]])*invAVcov[[j,k]]*(corrAVbootavg[[k,2]]-f[a,b,corrAVbootavg[[k,1]]]), {j, IPav, FPav}], {k, IPav, FPav}], {a,b}, Method -> "NelderMead"];
dof=(FPav-IPav+1)-2;
Abs[a]/. tmp[[2,1]], Abs[b] /. tmp[[2,2]],tmp[[1]]/dof};
(*Outputs {|a|,|b|,\[Chi]^2/dof}*)


(* ::Input:: *)
AVcovariance=Table[Table[Mean[Table[(corrAVboot[[j,i]]-corrAVbootavg[[j,2]])(corrAVboot[[k,i]]-corrAVbootavg[[k,2]]),{i,Nboot}]],{j,NT/2+1}],{k,NT/2+1}];
invAVcov=Inverse[AVcovariance];


(* ::Input:: *)
(*Find the longest sequence of t values for which the effective mass is both real and positive.*)
tmpMasses={};
For[i=1,i<=Dimensions[meffAV][[1]],++i,If[Element[meffAV[[i,2]],Reals],AppendTo[tmpMasses,meffAV[[i,;;]]]]]
setsplit=Split[tmpMasses[[;;,1]],#2-#1==1&];
lengths=Table[Length[setsplit[[i,;;]]],{i,Dimensions[setsplit][[1]]}];
place=Position[lengths,Max[lengths]][[1,1]];
window=setsplit[[place]]
(*If there are fewer than a specified number of points that satisfy the above, the fitting must be done by hand.*)
If[Length[window]<=0.4*Length[meffAV],{Print["Plateau fitting must be done manually!"],AutomaticAV=False;}]


(* ::Input:: *)
If[AutomaticAV&&UseAV,TInitial=window[[1]];TFinal=window[[-1]];
tmp=CurveFitAV[TInitial,TFinal];
chisq=tmp[[-1]];
For[t=TInitial+1,t<TFinal,++t,tmp=CurveFitAV[t,TFinal];dofav=(TFinal-t+1)-2;
If[tmp[[-1]]<chisq&&tmp[[-1]]>ChiSqCrit[dofav,0.90],chisq=tmp[[-1]],{TInitial=t-1,Break[]}]];TInitial=t-1;]


(* ::Input:: *)
If[AutomaticAV&&UseAV,{For[cert=0.90,cert<=0.95,cert+=0.05,{If[AutomaticAV&&UseAV,Variation=2;Tests={};
For[y=TFinal,y>=TFinal-Variation,--y,For[x=Max[1,TInitial-Variation],x<=y-MinimumPlateau,++x,
tmp=CurveFitAV[x,y];
dofav=(y-x+1)-2;
chisq=tmp[[-1]];
If[chisq>ChiSqCrit[dofav,cert],Continue[]];
chisq/=(y-x+1)^2;
AppendTo[Tests,{chisq,x,y}]]];
If[Length[Tests]!=0,{FitRange=Tests[[Position[Tests[[;;,1]],Min[Tests[[;;,1]]]][[1,1]]]],Break[]},Print["No valid fit ranges at cert = "<>ToString[cert]]];]}]}];
(*Select the fit range for which \[Chi]^2/dof(tf-ti+1)^2 is minimised. This attempts to simultaneously maximise the legnth of the fit range and minimise \[Chi]^2/dof. In some cases the initial certainty may be too restrictive. If so, the certainty is changed to widen the search.*)


(* ::Input:: *)
If[Length[Tests]==0,AutomaticAV=False];(*If cert is still too restrictive, the results have to be done by hand.*)


(* ::Input:: *)
If[AutomaticAV&&UseAV,{{IPav=FitRange[[2]],FPav=FitRange[[3]]},dofav=(FPav-IPav+1)-2}]
If[Not[AutomaticAV]&&UseAV,{{IPav=_SED_IPAV_,FPav=_SED_FPAV_},dofav=(FPav-IPav+1)-2}]


(* ::Input:: *)
ProgressIndicator[Dynamic[w3],{1,Nboot}]


(* ::Input:: *)
If[UseAV,Provisional=CurveFitAV[IPav,FPav]](*These are the provisional results that are subject to change during the robustness check.*)
(*Squared decay-constant, mass, \[Chi]^2/dof*)


(* ::Input:: *)
If[UseAV,BootsProvisional=BootstrapAV[IPav,FPav]];(*BootsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[UseAV,ErrsProvisional=BootErrorsAV[BootsProvisional]](*ErrsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[AutomaticAV&&UseAV,For[t=IPav+1,t<FPav+2-MinimumPlateau,++t,{temp=CurveFitAV[t,FPav];tempBoots=BootstrapAV[t,FPav];tempErr=BootErrorsAV[tempBoots];If[Abs[temp[[2]]-Provisional[[2]]]<ErrsProvisional[[2]],{AVFinal=Provisional,AVBoots=BootsProvisional,AVErr=ErrsProvisional,Break[]},{Provisional=temp,BootsProvisional=tempBoots,ErrsProvisional=tempErr}]}]]


(* ::Input:: *)
If[Not[AutomaticAV]&&UseAV,{AVFinal=Provisional,AVBoots=BootsProvisional,AVErr=ErrsProvisional}];


(* ::Input:: *)
{{IPav,FPav},dofav}


(* ::Input:: *)
If[UseAV,FINALSAV=ArrayReshape[Riffle[AVFinal,AVErr],{Length[AVFinal],2}]]
If[UseAV,FinalChiSq=FINALSAV[[-1,1]]];


(* ::Input:: *)
If[UseAV,For[i=1,i<=Nboot,++i,Write[StringJoin[Outputpath, "/", FilenameSlug,"_axialvector_masses_boots.csv"],(AVBoots[[i,2]]*GradientFlow)^2]], Write[StringJoin[Outputpath, "/", FilenameSlug,"_axialvector_masses_boots.csv"], ""]];
If[UseAV,For[i=1,i<=Nboot,++i,Write[StringJoin[Outputpath, "/", FilenameSlug,"_axialvector_decayconsts_boots.csv"],AVBoots[[i,1]]*(ZA*GradientFlow)^2]], Write[StringJoin[Outputpath, "/", FilenameSlug,"_axialvector_decayconsts_boots.csv"], ""]];
(*Write all bootstraps to an output file for use in the chiral extrapolation.*)


(* ::Input:: *)
If[UseAV,i=1;data=Abs[AVBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[AVFinal[[i]],AVErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[AVFinal[[i]]],", ",ToString[AVErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*SubsuperscriptBox[\(f\), \(AV\), \(2\)]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]
If[UseAV,i=2;data=Abs[AVBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[AVFinal[[i]],AVErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[AVFinal[[i]]],", ",ToString[AVErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*SubscriptBox[\(m\), \(AV\)]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]


(* ::Input:: *)
If[UseAV,FINALPLOTAV=Plot[{FINALSAV[[2,1]]+FINALSAV[[2,2]],FINALSAV[[2,1]]-FINALSAV[[2,2]]},{t,IPav-0.5,FPav+0.5},Filling->True];
If[AutomaticAV,FullPlotAV=Show[meffAVplot,FINALPLOTAV,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]]]]]]
If[Not[AutomaticAV],FullPlotAv=Show[meffAVplot,FINALPLOTAV,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]],"\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"Manual\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\")\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)"]]]


(* ::Input:: *)
If[UseAV,Export[StringJoin[Outputpath,"/axialvector.pdf"],FullPlotAV,"PDF"], Export[StringJoin[Outputpath,"/axialvector.pdf"], ""]];
If[UseAV,HistogramPlotAV=Histogram[AVBoots[[;;,-1]],{BinWidth},"Probability",ImageSize->600,ChartStyle->Darker[Blue],FrameLabel->{"\!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.","Normalised Frequency"},LabelStyle->Directive[20,Black],Frame->True,PlotRange->Automatic,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT], ", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=",ToString[NumberForm[m,{Digitsm,DPm}]]]]]


(* ::Input:: *)
If[UseAV,Export[StringJoin[Outputpath,"/axialvectorCSD.pdf"],HistogramPlotAV,"PDF"], Export[StringJoin[Outputpath,"/axialvectorCSD.pdf"], ""]];


(* ::Input:: *)
(*Square measurements, if necessary, accounting for propagation of errors and renormalise decay constants.*)
If[UseAV,fAV2=Around[FINALSAV[[1,1]],FINALSAV[[1,2]]];
massAV=Around[FINALSAV[[2,1]],FINALSAV[[2,2]]];
fAV2=fAV2*(GradientFlow)^2;
massAV=(GradientFlow*massAV)^2;
FINALSAV[[1,1]]=ZA^2*fAV2["Value"];
FINALSAV[[1,2]]=ZA^2*fAV2["Error"];
FINALSAV[[2,1]]=massAV["Value"];
FINALSAV[[2,2]]=massAV["Error"];]


(* ::Input:: *)
If[UseAV,Export[StringJoin[Outputpath,"/output_axialvector.txt"],{Flatten[FINALSAV],ToString[IPav],ToString[FPav]},"Table","FieldSeparators"->"\n"], Export[StringJoin[Outputpath,"/output_axialvector.txt"], ""]];


(* ::Section:: *)
(* Scalar *)


(* ::Input:: *)
UseS=_SED_USES_;
AutomaticS=_SED_AUTOMATICS_;
MinimumPlateau=3;


(* ::Input:: *)
ch=2;
Idimport=Table[RAWDATA[[Length[channels]*(i-1)+ch]],{i,1,Cnfg}];


(* ::Input:: *)
For[i=1,i<=Cnfg,i++,Idimport[[i]]=Append[Idimport[[i]],Idimport[[i,3]]]];


(* ::Input:: *)
corrS=Table[Table[NS^3(Idimport[[j,i+2]]+Idimport[[j,NT+4-i]])/2,{i,1,NT/2+1}],{j,1,Cnfg}];


(* ::Input:: *)
corrS=Transpose[corrS];


(* ::Input:: *)
RandInt=Table[Table[Round[RandomReal[{0,1}]*(Cnfg-1)+1],{i,Nboot}],{j,Cnfg}];
corrSboot=Table[Table[Mean[Table[corrS[[k,RandInt[[i,j]]]],{i,Cnfg}]],{j,Nboot}],{k,NT/2+1}];
corrSbootavg=Table[{k-1,Mean[corrSboot[[k]]],StandardDeviation[corrSboot[[k]]]},{k,NT/2+1}];


(* ::Input:: *)
meffS=Table[{i,Mean[Table[ArcCosh[(corrSboot[[i,j]]+corrSboot[[i+2,j]])/(2*corrSboot[[i+1,j]])],{j,Nboot}]],StandardDeviation[Table[ArcCosh[(corrSboot[[i,j]]+corrSboot[[i+2,j]])/(2*corrSboot[[i+1,j]])],{j,Nboot}]]},{i,1,NT/2-1}];


(* ::Input:: *)
meffSplot=ErrorListPlot[Table[{{meffS[[t,1]],meffS[[t,2]]},ErrorBar[meffS[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->Automatic, FrameLabel->{"t","\!\(\*SubsuperscriptBox[\(am\), \(eff\), \(S\)]\)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->Automatic,AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[0,0,1],PlotLabel->StringJoin["\[DoubleStruckOne]. Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]]]]


(* ::Input:: *)
Export[StringJoin[Outputpath,"/scalar.pdf"],meffSplot,"PDF"];


(* ::Input:: *)
Scorrplot=ErrorListPlot[Table[{{corrSbootavg[[t,1]],corrSbootavg[[t,2]]},ErrorBar[corrSbootavg[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->{"o",15}, FrameLabel->{"t","\!\(\*SubscriptBox[\(C\), \(\[DoubleStruckOne]\)]\)(t)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->{{0,NT},All},AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[1,0,0]]


(* ::Input:: *)
Clear[f]
f[a_,b_,t_]:=(a/b)*(Exp[-b*t]+Exp[-b(NT-t)]);
(*Fitting function for the matrix element (a) and the mass (b)*)
BootstrapS[IPs_,FPs_]:=Table[w4=i;tmp = NMinimize[Sum[Sum[(corrSboot[[j,i]]-f[a,b,corrSbootavg[[j,1]]])*invScov[[j,k]]*(corrSboot[[k,i]]-f[a,b,corrSbootavg[[k,1]]]), {j, IPs, FPs}], {k, IPs, FPs}], {a,b}, Method -> "NelderMead"];
dof=(FPs-IPs+1)-2;
 {a/. tmp[[2,1]], b /. tmp[[2,2]],tmp[[1]]/dof}, {i, 1, Nboot}];
(*Outputs {{a,b,\[Chi]^2/dof},...} Nboot times*)
BootErrorsS[input_]:=Table[DistributionFitTest[Abs[input[[;;,i]]],Automatic,"HypothesisTestData"]["FittedDistribution"][[2]],{i,1,Dimensions[input][[2]]}]
CurveFitS[IPs_,FPs_]:={tmp = NMinimize[Sum[Sum[(corrSbootavg[[j,2]]-f[a,b,corrSbootavg[[j,1]]])*invScov[[j,k]]*(corrSbootavg[[k,2]]-f[a,b,corrSbootavg[[k,1]]]), {j, IPs, FPs}], {k, IPs, FPs}], {a,b}, Method -> "NelderMead"];
dof=(FPs-IPs+1)-2;
Abs[a]/. tmp[[2,1]], Abs[b] /. tmp[[2,2]],tmp[[1]]/dof};
(*Outputs {|a|,|b|,\[Chi]^2/dof}*)


(* ::Input:: *)
Scovariance=Table[Table[Mean[Table[(corrSboot[[j,i]]-corrSbootavg[[j,2]])(corrSboot[[k,i]]-corrSbootavg[[k,2]]),{i,Nboot}]],{j,NT/2+1}],{k,NT/2+1}];
invScov=Inverse[Scovariance];


(* ::Input:: *)
(*Find the longest sequence of t values for which the effective mass is both real and positive.*)
tmpMasses={};
For[i=1,i<=Dimensions[meffS][[1]],++i,If[Element[meffS[[i,2]],Reals],AppendTo[tmpMasses,meffS[[i,;;]]]]]
setsplit=Split[tmpMasses[[;;,1]],#2-#1==1&];
lengths=Table[Length[setsplit[[i,;;]]],{i,Dimensions[setsplit][[1]]}];
place=Position[lengths,Max[lengths]][[1,1]];
window=setsplit[[place]]
(*If there are fewer than a specified number of points that satisfy the above, the fitting must be done by hand.*)
If[Length[window]<=0.4*Length[meffS],{Print["Plateau fitting must be done manually!"],AutomaticS=False;}]


(* ::Input:: *)
If[AutomaticS&&UseS,TInitial=window[[1]];TFinal=window[[-1]];
tmp=CurveFitS[TInitial,TFinal];
chisq=tmp[[-1]];
For[t=TInitial+1,t<TFinal,++t,tmp=CurveFitS[t,TFinal];dofs=(TFinal-t+1)-2;
If[tmp[[-1]]<chisq&&tmp[[-1]]>ChiSqCrit[dofs,0.90],chisq=tmp[[-1]],{TInitial=t-1,Break[]}]];TInitial=t-1;]


(* ::Input:: *)
If[AutomaticS&&UseS,{For[cert=0.90,cert<=0.95,cert+=0.05,{If[AutomaticS&&UseS,Variation=2;Tests={};
For[y=TFinal,y>=TFinal-Variation,--y,For[x=Max[1,TInitial-Variation],x<=y-MinimumPlateau,++x,
tmp=CurveFitS[x,y];
dofs=(y-x+1)-2;
chisq=tmp[[-1]];
If[chisq>ChiSqCrit[dofs,cert],Continue[]];
chisq/=(y-x+1)^2;
AppendTo[Tests,{chisq,x,y}]]];
If[Length[Tests]!=0,{FitRange=Tests[[Position[Tests[[;;,1]],Min[Tests[[;;,1]]]][[1,1]]]],Break[]},Print["No valid fit ranges at cert = "<>ToString[cert]]];]}]}];
(*Select the fit range for which \[Chi]^2/dof(tf-ti+1)^2 is minimised. This attempts to simultaneously maximise the legnth of the fit range and minimise \[Chi]^2/dof. In some cases the initial certainty may be too restrictive. If so, the certainty is changed to widen the search.*)


(* ::Input:: *)
If[Length[Tests]==0,AutomaticS=False];(*If cert is still too restrictive, the results have to be done by hand.*)


(* ::Input:: *)
If[AutomaticS&&UseS,{{IPs=FitRange[[2]],FPs=FitRange[[3]]},dofs=(FPs-IPs+1)-2}]
If[Not[AutomaticS]&&UseS,{{IPs=_SED_IPS_,FPs=_SED_FPS_},dofs=(FPs-IPs+1)-2}]


(* ::Input:: *)
ProgressIndicator[Dynamic[w4],{1,Nboot}]


(* ::Input:: *)
If[UseS,Provisional=CurveFitS[IPs,FPs]](*These are the provisional results that are subject to change during the robustness check.*)
(*Squared matrix-element, mass, \[Chi]^2/dof*)


(* ::Input:: *)
If[UseS,BootsProvisional=BootstrapS[IPs,FPs]];(*BootsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[UseS,ErrsProvisional=BootErrorsS[BootsProvisional]](*ErrsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[AutomaticS&&UseS,For[t=IPs+1,t<FPs+2-MinimumPlateau,++t,{temp=CurveFitS[t,FPs];tempBoots=BootstrapS[t,FPs];tempErr=BootErrorsS[tempBoots];If[Abs[temp[[2]]-Provisional[[2]]]<ErrsProvisional[[2]],{SFinal=Provisional,SBoots=BootsProvisional,SErr=ErrsProvisional,Break[]},{Provisional=temp,BootsProvisional=tempBoots,ErrsProvisional=tempErr}]}]]


(* ::Input:: *)
If[Not[AutomaticS]&&UseS,{SFinal=Provisional,SBoots=BootsProvisional,SErr=ErrsProvisional}];


(* ::Input:: *)
{{IPs,FPs},dofs}


(* ::Input:: *)
If[UseS,FINALSS=ArrayReshape[Riffle[SFinal,SErr],{Length[SFinal],2}]]
If[UseS,FinalChiSq=FINALSS[[-1,1]]];


(* ::Input:: *)
If[UseS,For[i=1,i<=Nboot,++i,Write[StringJoin[Outputpath, "/", FilenameSlug,"_scalar_masses_boots.csv"],(SBoots[[i,2]]*GradientFlow)^2]], Write[StringJoin[Outputpath, "/", FilenameSlug,"_scalar_masses_boots.csv"], ""]];
(*Write all bootstraps to an output file for use in the chiral extrapolation.*)


(* ::Input:: *)
If[UseS,i=1;data=Abs[SBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[SFinal[[i]],SErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[SFinal[[i]]],", ",ToString[SErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*TemplateBox[{\"0\"},\n\"Bra\"]\)\!\(\*SubscriptBox[\(\[ScriptCapitalO]\), \(S\)]\)\!\(\*TemplateBox[{\"S\"},\n\"Ket\"]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\[DoubleStruckOne]. Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]
If[UseS,i=2;data=Abs[SBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[SFinal[[i]],SErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[SFinal[[i]]],", ",ToString[SErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*SubscriptBox[\(m\), \(S\)]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\[DoubleStruckOne]. Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]


(* ::Input:: *)
If[UseS,FINALPLOTS=Plot[{FINALSS[[2,1]]+FINALSS[[2,2]],FINALSS[[2,1]]-FINALSS[[2,2]]},{t,IPs-0.5,FPs+0.5},Filling->True];
If[AutomaticS,FullPlotS=Show[meffSplot,FINALPLOTS,PlotLabel->StringJoin["\[DoubleStruckOne]. Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]]]]]]
If[Not[AutomaticS],FullPlotS=Show[meffSplot,FINALPLOTS,PlotLabel->StringJoin["\[DoubleStruckOne]. Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]],"\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"Manual\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\")\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)"]]]


(* ::Input:: *)
If[UseS,Export[StringJoin[Outputpath,"/scalar.pdf"],FullPlotS,"PDF"], Export[StringJoin[Outputpath,"/scalar.pdf"], ""]];
If[UseS,HistogramPlotS=Histogram[SBoots[[;;,-1]],{BinWidth},"Probability",ImageSize->600,ChartStyle->Darker[Blue],FrameLabel->{"\!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.","Normalised Frequency"},LabelStyle->Directive[20,Black],Frame->True,PlotRange->Automatic,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\[DoubleStruckOne]. Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT], ", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=",ToString[NumberForm[m,{Digitsm,DPm}]]]]]


(* ::Input:: *)
If[UseS,Export[StringJoin[Outputpath,"/scalarCSD.pdf"],HistogramPlotS,"PDF"], Export[StringJoin[Outputpath,"/scalarCSD.pdf"], ""]];


(* ::Input:: *)
(*Square measurements, if necessary, accounting for propagation of errors and renormalise decay constants.*)
If[UseS,massS=Around[FINALSS[[2,1]],FINALSS[[2,2]]];
massS=(GradientFlow*massS)^2;
FINALSS[[2,1]]=massS["Value"];
FINALSS[[2,2]]=massS["Error"];]


(* ::Input:: *)
If[UseS,Export[StringJoin[Outputpath,"/output_scalar.txt"],{Flatten[FINALSS],ToString[IPs],ToString[FPs]},"Table","FieldSeparators"->"\n"], Export[StringJoin[Outputpath,"/output_scalar.txt"], ""]];


(* ::Section:: *)
(* Tensor *)


(* ::Input:: *)
UseT=_SED_USET_;
AutomaticT=_SED_AUTOMATICT_;
MinimumPlateau=3;


(* ::Input:: *)
ch=7;
g0g1import=Table[RAWDATA[[Length[channels]*(i-1)+ch]],{i,1,Cnfg}];
g0g2import=Table[RAWDATA[[Length[channels]*(i-1)+ch+1]],{i,1,Cnfg}];
g0g3import=Table[RAWDATA[[Length[channels]*(i-1)+ch+2]],{i,1,Cnfg}];


(* ::Input:: *)
For[i=1,i<=Cnfg,i++,g0g1import[[i]]=Append[g0g1import[[i]],g0g1import[[i,3]]]];
For[i=1,i<=Cnfg,i++,g0g2import[[i]]=Append[g0g2import[[i]],g0g2import[[i,3]]]];
For[i=1,i<=Cnfg,i++,g0g3import[[i]]=Append[g0g3import[[i]],g0g3import[[i,3]]]];


(* ::Input:: *)
corrT=Table[Table[NS^3(g0g1import[[j,i+2]]+g0g2import[[j,i+2]]+g0g2import[[j,i+2]]+g0g1import[[j,NT+4-i]]+g0g2import[[j,NT+4-i]]+g0g3import[[j,NT+4-i]])/6,{i,1,NT/2+1}],{j,1,Cnfg}];


(* ::Input:: *)
corrT=Transpose[corrT];


(* ::Input:: *)
RandInt=Table[Table[Round[RandomReal[{0,1}]*(Cnfg-1)+1],{i,Nboot}],{j,Cnfg}];
corrTboot=Table[Table[Mean[Table[corrT[[k,RandInt[[i,j]]]],{i,Cnfg}]],{j,Nboot}],{k,NT/2+1}];
corrTbootavg=Table[{k-1,Mean[corrTboot[[k]]],StandardDeviation[corrTboot[[k]]]},{k,NT/2+1}];


(* ::Input:: *)
meffT=Table[{i,Mean[Table[ArcCosh[(corrTboot[[i,j]]+corrTboot[[i+2,j]])/(2*corrTboot[[i+1,j]])],{j,Nboot}]],StandardDeviation[Table[ArcCosh[(corrTboot[[i,j]]+corrTboot[[i+2,j]])/(2*corrTboot[[i+1,j]])],{j,Nboot}]]},{i,1,NT/2-1}];


(* ::Input:: *)
meffTplot=ErrorListPlot[Table[{{meffT[[t,1]],meffT[[t,2]]},ErrorBar[meffT[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->Automatic, FrameLabel->{"t","\!\(\*SubsuperscriptBox[\(am\), \(eff\), \(T\)]\)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->Automatic,AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[0,0,1],PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]]]]


(* ::Input:: *)
Export[StringJoin[Outputpath,"/tensor.pdf"],meffTplot,"PDF"];


(* ::Input:: *)
Tcorrplot=ErrorListPlot[Table[{{corrTbootavg[[t,1]],corrTbootavg[[t,2]]},ErrorBar[corrTbootavg[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->{"o",15}, FrameLabel->{"t","\!\(\*SubscriptBox[\(C\), \(\[Gamma]0\[Gamma]k\)]\)(t)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->{{0,NT},All},AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[1,0,0]]


(* ::Input:: *)
Clear[f]
f[a_,b_,t_]:=(a/b)*(Exp[-b*t]+Exp[-b(NT-t)]);
(*Fitting function for the matrix element (a) and the mass (b)*)
BootstrapT[IPt_,FPt_]:=Table[w5=i;tmp = NMinimize[Sum[Sum[(corrTboot[[j,i]]-f[a,b,corrTbootavg[[j,1]]])*invTcov[[j,k]]*(corrTboot[[k,i]]-f[a,b,corrTbootavg[[k,1]]]), {j, IPt, FPt}], {k, IPt, FPt}], {a,b}, Method -> "NelderMead"];
dof=(FPt-IPt+1)-2;
 {a/. tmp[[2,1]], b /. tmp[[2,2]],tmp[[1]]/dof}, {i, 1, Nboot}];
(*Outputs {{a,b,\[Chi]^2/dof},...} Nboot times*)
BootErrorsT[input_]:=Table[DistributionFitTest[Abs[input[[;;,i]]],Automatic,"HypothesisTestData"]["FittedDistribution"][[2]],{i,1,Dimensions[input][[2]]}]
CurveFitT[IPt_,FPt_]:={tmp = NMinimize[Sum[Sum[(corrTbootavg[[j,2]]-f[a,b,corrTbootavg[[j,1]]])*invTcov[[j,k]]*(corrTbootavg[[k,2]]-f[a,b,corrTbootavg[[k,1]]]), {j, IPt, FPt}], {k, IPt, FPt}], {a,b}, Method -> "NelderMead"];
dof=(FPt-IPt+1)-2;
Abs[a]/. tmp[[2,1]], Abs[b] /. tmp[[2,2]],tmp[[1]]/dof};
(*Outputs {|a|,|b|,\[Chi]^2/dof}*)


(* ::Input:: *)
Tcovariance=Table[Table[Mean[Table[(corrTboot[[j,i]]-corrTbootavg[[j,2]])(corrTboot[[k,i]]-corrTbootavg[[k,2]]),{i,Nboot}]],{j,NT/2+1}],{k,NT/2+1}];
invTcov=Inverse[Tcovariance];


(* ::Input:: *)
(*Find the longest sequence of t values for which the effective mass is both real and positive.*)
tmpMasses={};
For[i=1,i<=Dimensions[meffT][[1]],++i,If[Element[meffT[[i,2]],Reals],AppendTo[tmpMasses,meffT[[i,;;]]]]]
setsplit=Split[tmpMasses[[;;,1]],#2-#1==1&];
lengths=Table[Length[setsplit[[i,;;]]],{i,Dimensions[setsplit][[1]]}];
place=Position[lengths,Max[lengths]][[1,1]];
window=setsplit[[place]]
(*If there are fewer than a specified number of points that satisfy the above, the fitting must be done by hand.*)
If[Length[window]<=0.4*Length[meffT],{Print["Plateau fitting must be done manually!"],AutomaticT=False;}]


(* ::Input:: *)
If[AutomaticT&&UseT,TInitial=window[[1]];TFinal=window[[-1]];
tmp=CurveFitT[TInitial,TFinal];
chisq=tmp[[-1]];
For[t=TInitial+1,t<TFinal,++t,tmp=CurveFitT[t,TFinal];doft=(TFinal-t+1)-2;
If[tmp[[-1]]<chisq&&tmp[[-1]]>ChiSqCrit[doft,0.90],chisq=tmp[[-1]],{TInitial=t-1,Break[]}]];TInitial=t-1;]


(* ::Input:: *)
If[AutomaticT&&UseT,{For[cert=0.90,cert<=0.95,cert+=0.05,{If[AutomaticT&&UseT,Variation=2;Tests={};
For[y=TFinal,y>=TFinal-Variation,--y,For[x=Max[1,TInitial-Variation],x<=y-MinimumPlateau,++x,
tmp=CurveFitT[x,y];
doft=(y-x+1)-2;
chisq=tmp[[-1]];
If[chisq>ChiSqCrit[doft,cert],Continue[]];
chisq/=(y-x+1)^2;
AppendTo[Tests,{chisq,x,y}]]];
If[Length[Tests]!=0,{FitRange=Tests[[Position[Tests[[;;,1]],Min[Tests[[;;,1]]]][[1,1]]]],Break[]},Print["No valid fit ranges at cert = "<>ToString[cert]]];]}]}];
(*Select the fit range for which \[Chi]^2/dof(tf-ti+1)^2 is minimised. This attempts to simultaneously maximise the legnth of the fit range and minimise \[Chi]^2/dof. In some cases the initial certainty may be too restrictive. If so, the certainty is changed to widen the search.*)


(* ::Input:: *)
If[Length[Tests]==0,AutomaticT=False];(*If cert is still too restrictive, the results have to be done by hand.*)


(* ::Input:: *)
If[AutomaticT&&UseT,{{IPt=FitRange[[2]],FPt=FitRange[[3]]},doft=(FPt-IPt+1)-2}]
If[Not[AutomaticT]&&UseT,{{IPt=_SED_IPT_,FPt=_SED_FPT_},doft=(FPt-IPt+1)-2}]


(* ::Input:: *)
ProgressIndicator[Dynamic[w5],{1,Nboot}]


(* ::Input:: *)
If[UseT,Provisional=CurveFitT[IPt,FPt]](*These are the provisional results that are subject to change during the robustness check.*)
(*Squared matrix-element, mass, \[Chi]^2/dof*)


(* ::Input:: *)
If[UseT,BootsProvisional=BootstrapT[IPt,FPt]];(*BootsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[UseT,ErrsProvisional=BootErrorsT[BootsProvisional]](*ErrsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[AutomaticT&&UseT,For[t=IPt+1,t<FPt+2-MinimumPlateau,++t,{temp=CurveFitT[t,FPt];tempBoots=BootstrapT[t,FPt];tempErr=BootErrorsT[tempBoots];If[Abs[temp[[2]]-Provisional[[2]]]<ErrsProvisional[[2]],{TenFinal=Provisional,TBoots=BootsProvisional,TErr=ErrsProvisional,Break[]},{Provisional=temp,BootsProvisional=tempBoots,ErrsProvisional=tempErr}]}]]


(* ::Input:: *)
If[Not[AutomaticT]&&UseT,{TenFinal=Provisional,TBoots=BootsProvisional,TErr=ErrsProvisional}];


(* ::Input:: *)
{{IPt,FPt},doft}


(* ::Input:: *)
If[UseT,FINALST=ArrayReshape[Riffle[TenFinal,TErr],{Length[TenFinal],2}]]
If[UseT,FinalChiSq=FINALST[[-1,1]]];


(* ::Input:: *)
If[UseT,For[i=1,i<=Nboot,++i,Write[StringJoin[Outputpath, "/", FilenameSlug,"_tensor_masses_boots.csv"],(TBoots[[i,2]]*GradientFlow)^2]], Write[StringJoin[Outputpath, "/", FilenameSlug,"_tensor_masses_boots.csv"], ""]];
(*Write all bootstraps to an output file for use in the chiral extrapolation.*)


(* ::Input:: *)
If[UseT,i=1;data=Abs[TBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[TenFinal[[i]],TErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[TenFinal[[i]]],", ",ToString[TErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*TemplateBox[{\"0\"},\n\"Bra\"]\)\!\(\*SubscriptBox[\(\[ScriptCapitalO]\), \(T\)]\)\!\(\*TemplateBox[{\"S\"},\n\"Ket\"]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]
If[UseT,i=2;data=Abs[TBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[TenFinal[[i]],TErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[TenFinal[[i]]],", ",ToString[TErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*SubscriptBox[\(m\), \(T\)]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]


(* ::Input:: *)
If[UseT,FINALPLOTT=Plot[{FINALST[[2,1]]+FINALST[[2,2]],FINALST[[2,1]]-FINALST[[2,2]]},{t,IPt-0.5,FPt+0.5},Filling->True];
If[AutomaticT,FullPlotT=Show[meffTplot,FINALPLOTT,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]]]]]]
If[Not[AutomaticT],FullPlotT=Show[meffTplot,FINALPLOTT,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]],"\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"Manual\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\")\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)"]]]


(* ::Input:: *)
If[UseT,Export[StringJoin[Outputpath,"/tensor.pdf"],FullPlotT,"PDF"], Export[StringJoin[Outputpath,"/tensor.pdf"], ""]];
If[UseT,HistogramPlotT=Histogram[TBoots[[;;,-1]],{BinWidth},"Probability",ImageSize->600,ChartStyle->Darker[Blue],FrameLabel->{"\!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.","Normalised Frequency"},LabelStyle->Directive[20,Black],Frame->True,PlotRange->Automatic,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT], ", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=",ToString[NumberForm[m,{Digitsm,DPm}]]]]]


(* ::Input:: *)
If[UseT,Export[StringJoin[Outputpath,"/tensorCSD.pdf"],HistogramPlotT,"PDF"], Export[StringJoin[Outputpath,"/tensorCSD.pdf"], ""]];


(* ::Input:: *)
(*Square measurements, if necessary, accounting for propagation of errors and renormalise decay constants.*)
If[UseT,massT=Around[FINALST[[2,1]],FINALST[[2,2]]];
massT=(GradientFlow*massT)^2;
FINALST[[2,1]]=massT["Value"];
FINALST[[2,2]]=massT["Error"];]


(* ::Input:: *)
If[UseT,Export[StringJoin[Outputpath,"/output_tensor.txt"],{Flatten[FINALST],ToString[IPt],ToString[FPt]},"Table","FieldSeparators"->"\n"], Export[StringJoin[Outputpath,"/output_tensor.txt"], ""]];


(* ::Section:: *)
(* Axial-Tensor *)


(* ::Input:: *)
UseAT=_SED_USEAT_;
AutomaticAT=_SED_AUTOMATICAT_;
MinimumPlateau=3;


(* ::Input:: *)
ch=14;
g5g0g1import=Table[RAWDATA[[Length[channels]*(i-1)+ch]],{i,1,Cnfg}];
g5g0g2import=Table[RAWDATA[[Length[channels]*(i-1)+ch+1]],{i,1,Cnfg}];
g5g0g3import=Table[RAWDATA[[Length[channels]*(i-1)+ch+2]],{i,1,Cnfg}];


(* ::Input:: *)
For[i=1,i<=Cnfg,i++,g5g0g1import[[i]]=Append[g5g0g1import[[i]],g5g0g1import[[i,3]]]];
For[i=1,i<=Cnfg,i++,g5g0g2import[[i]]=Append[g5g0g2import[[i]],g5g0g2import[[i,3]]]];
For[i=1,i<=Cnfg,i++,g5g0g3import[[i]]=Append[g5g0g3import[[i]],g5g0g3import[[i,3]]]];


(* ::Input:: *)
corrAT=Table[Table[NS^3(g5g0g1import[[j,i+2]]+g5g0g2import[[j,i+2]]+g5g0g2import[[j,i+2]]+g5g0g1import[[j,NT+4-i]]+g5g0g2import[[j,NT+4-i]]+g5g0g3import[[j,NT+4-i]])/6,{i,1,NT/2+1}],{j,1,Cnfg}];


(* ::Input:: *)
corrAT=Transpose[corrAT];


(* ::Input:: *)
RandInt=Table[Table[Round[RandomReal[{0,1}]*(Cnfg-1)+1],{i,Nboot}],{j,Cnfg}];
corrATboot=Table[Table[Mean[Table[corrAT[[k,RandInt[[i,j]]]],{i,Cnfg}]],{j,Nboot}],{k,NT/2+1}];
corrATbootavg=Table[{k-1,Mean[corrATboot[[k]]],StandardDeviation[corrATboot[[k]]]},{k,NT/2+1}];


(* ::Input:: *)
meffAT=Table[{i,Mean[Table[ArcCosh[(corrATboot[[i,j]]+corrATboot[[i+2,j]])/(2*corrATboot[[i+1,j]])],{j,Nboot}]],StandardDeviation[Table[ArcCosh[(corrATboot[[i,j]]+corrATboot[[i+2,j]])/(2*corrATboot[[i+1,j]])],{j,Nboot}]]},{i,1,NT/2-1}];


(* ::Input:: *)
meffATplot=ErrorListPlot[Table[{{meffAT[[t,1]],meffAT[[t,2]]},ErrorBar[meffAT[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->Automatic, FrameLabel->{"t","\!\(\*SubsuperscriptBox[\(am\), \(eff\), \(AT\)]\)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->Automatic,AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[0,0,1],PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]]]]


(* ::Input:: *)
Export[StringJoin[Outputpath,"/atensor.pdf"],meffATplot,"PDF"];


(* ::Input:: *)
ATcorrplot=ErrorListPlot[Table[{{corrATbootavg[[t,1]],corrATbootavg[[t,2]]},ErrorBar[corrATbootavg[[t,3]]]},{t,1,NT/2-1}], ImageSize->600, PlotMarkers->{"o",15}, FrameLabel->{"t","\!\(\*SubscriptBox[\(C\), \(\[Gamma]5\[Gamma]0\[Gamma]k\)]\)(t)"},LabelStyle->Directive[15,Black],Frame->True,PlotRange->{{0,NT},All},AxesStyle->Thick,AxesStyle->Black,PlotStyle->RGBColor[1,0,0]]


(* ::Input:: *)
Clear[f]
f[a_,b_,t_]:=(a/b)*(Exp[-b*t]+Exp[-b(NT-t)]);
(*Fitting function for the matrix element (a) and the mass (b)*)
BootstrapAT[IPat_,FPat_]:=Table[w6=i;tmp = NMinimize[Sum[Sum[(corrATboot[[j,i]]-f[a,b,corrATbootavg[[j,1]]])*invATcov[[j,k]]*(corrATboot[[k,i]]-f[a,b,corrATbootavg[[k,1]]]), {j, IPat, FPat}], {k, IPat, FPat}], {a,b}, Method -> "NelderMead"];
dof=(FPat-IPat+1)-2;
 {a/. tmp[[2,1]], b /. tmp[[2,2]],tmp[[1]]/dof}, {i, 1, Nboot}];
(*Outputs {{a,b,\[Chi]^2/dof},...} Nboot times*)
BootErrorsAT[input_]:=Table[DistributionFitTest[Abs[input[[;;,i]]],Automatic,"HypothesisTestData"]["FittedDistribution"][[2]],{i,1,Dimensions[input][[2]]}]
CurveFitAT[IPat_,FPat_]:={tmp = NMinimize[Sum[Sum[(corrATbootavg[[j,2]]-f[a,b,corrATbootavg[[j,1]]])*invATcov[[j,k]]*(corrATbootavg[[k,2]]-f[a,b,corrATbootavg[[k,1]]]), {j, IPat, FPat}], {k, IPat, FPat}], {a,b}, Method -> "NelderMead"];
dof=(FPat-IPat+1)-2;
Abs[a]/. tmp[[2,1]], Abs[b] /. tmp[[2,2]],tmp[[1]]/dof};
(*Outputs {|a|,|b|,\[Chi]^2/dof}*)


(* ::Input:: *)
ATcovariance=Table[Table[Mean[Table[(corrATboot[[j,i]]-corrATbootavg[[j,2]])(corrATboot[[k,i]]-corrATbootavg[[k,2]]),{i,Nboot}]],{j,NT/2+1}],{k,NT/2+1}];
invATcov=Inverse[ATcovariance];


(* ::Input:: *)
(*Find the longest sequence of t values for which the effective mass is both real and positive.*)
tmpMasses={};
For[i=1,i<=Dimensions[meffAT][[1]],++i,If[Element[meffAT[[i,2]],Reals],AppendTo[tmpMasses,meffAT[[i,;;]]]]]
setsplit=Split[tmpMasses[[;;,1]],#2-#1==1&];
lengths=Table[Length[setsplit[[i,;;]]],{i,Dimensions[setsplit][[1]]}];
place=Position[lengths,Max[lengths]][[1,1]];
window=setsplit[[place]]
(*If there are fewer than a specified number of points that satisfy the above, the fitting must be done by hand.*)
If[Length[window]<=0.4*Length[meffAT],{Print["Plateau fitting must be done manually!"],AutomaticAT=False;}]


(* ::Input:: *)
If[AutomaticAT&&UseAT,TInitial=window[[1]];TFinal=window[[-1]];
tmp=CurveFitAT[TInitial,TFinal];
chisq=tmp[[-1]];
For[t=TInitial+1,t<TFinal,++t,tmp=CurveFitAT[t,TFinal];dofat=(TFinal-t+1)-2;
If[tmp[[-1]]<chisq&&tmp[[-1]]>ChiSqCrit[dofat,0.90],chisq=tmp[[-1]],{TInitial=t-1,Break[]}]];TInitial=t-1;]


(* ::Input:: *)
If[AutomaticAT&&UseAT,{For[cert=0.90,cert<=0.95,cert+=0.05,{If[AutomaticAT&&UseAT,Variation=2;Tests={};
For[y=TFinal,y>=TFinal-Variation,--y,For[x=Max[1,TInitial-Variation],x<=y-MinimumPlateau,++x,
tmp=CurveFitAT[x,y];
dofat=(y-x+1)-2;
chisq=tmp[[-1]];
If[chisq>ChiSqCrit[dofat,cert],Continue[]];
chisq/=(y-x+1)^2;
AppendTo[Tests,{chisq,x,y}]]];
If[Length[Tests]!=0,{FitRange=Tests[[Position[Tests[[;;,1]],Min[Tests[[;;,1]]]][[1,1]]]],Break[]},Print["No valid fit ranges at cert = "<>ToString[cert]]];]}]}];
(*Select the fit range for which \[Chi]^2/dof(tf-ti+1)^2 is minimised. This attempts to simultaneously maximise the legnth of the fit range and minimise \[Chi]^2/dof. In some cases the initial certainty may be too restrictive. If so, the certainty is changed to widen the search.*)


(* ::Input:: *)
If[Length[Tests]==0,AutomaticAT=False];(*If cert is still too restrictive, the results have to be done by hand.*)


(* ::Input:: *)
If[AutomaticAT&&UseAT,{{IPat=FitRange[[2]],FPat=FitRange[[3]]},dofat=(FPat-IPat+1)-2}]
If[Not[AutomaticAT]&&UseAT,{{IPat=_SED_IPAT_,FPat=_SED_FPAT_},dofat=(FPat-IPat+1)-2}]


(* ::Input:: *)
ProgressIndicator[Dynamic[w6],{1,Nboot}]


(* ::Input:: *)
If[UseAT,Provisional=CurveFitAT[IPat,FPat]](*These are the provisional results that are subject to change during the robustness check.*)
(*Squared matrix-element, mass, \[Chi]^2/dof*)


(* ::Input:: *)
If[UseAT,BootsProvisional=BootstrapAT[IPat,FPat]];(*BootsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[UseAT,ErrsProvisional=BootErrorsAT[BootsProvisional]](*ErrsProvisional is used as a placeholder for the robustness check.*)


(* ::Input:: *)
If[AutomaticAT&&UseAT,For[t=IPat+1,t<FPat+2-MinimumPlateau,++t,{temp=CurveFitAT[t,FPat];tempBoots=BootstrapAT[t,FPat];tempErr=BootErrorsAT[tempBoots];If[Abs[temp[[2]]-Provisional[[2]]]<ErrsProvisional[[2]],{ATFinal=Provisional,ATBoots=BootsProvisional,ATErr=ErrsProvisional,Break[]},{Provisional=temp,BootsProvisional=tempBoots,ErrsProvisional=tempErr}]}]]


(* ::Input:: *)
If[Not[AutomaticAT]&&UseAT,{ATFinal=Provisional,ATBoots=BootsProvisional,ATErr=ErrsProvisional}];


(* ::Input:: *)
{{IPat,FPat},dofat}


(* ::Input:: *)
If[UseAT,FINALSAT=ArrayReshape[Riffle[ATFinal,ATErr],{Length[ATFinal],2}]]
If[UseAT,FinalChiSq=FINALSAT[[-1,1]]];


(* ::Input:: *)
If[UseAT,For[i=1,i<=Nboot,++i,Write[StringJoin[Outputpath, "/", FilenameSlug,"_axialtensor_masses_boots.csv"],(ATBoots[[i,2]]*GradientFlow)^2]], Write[StringJoin[Outputpath, "/", FilenameSlug,"_axialtensor_masses_boots.csv"], ""]];
(*Write all bootstraps to an output file for use in the chiral extrapolation.*)


(* ::Input:: *)
If[UseAT,i=1;data=Abs[ATBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[ATFinal[[i]],ATErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[ATFinal[[i]]],", ",ToString[ATErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*TemplateBox[{\"0\"},\n\"Bra\"]\)\!\(\*SubscriptBox[\(\[ScriptCapitalO]\), \(AT\)]\)\!\(\*TemplateBox[{\"S\"},\n\"Ket\"]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]
If[UseAT,i=2;data=Abs[ATBoots[[;;,i]]];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[NormalDistribution[ATFinal[[i]],ATErr[[i]]],x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[ATFinal[[i]]],", ",ToString[ATErr[[i]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{"\!\(\*SubscriptBox[\(m\), \(AT\)]\)","PDF"},LabelStyle->Directive[15,Black],Frame->True,ImageSize->600,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(", ToString[Nc], "), ", ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Cross]",ToString[NT],",  \[Beta]=", ToString[NumberForm[beta,{Digitsb,DPb}]], ", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm, DPm}]],". Bootstrap distribution"]]]


(* ::Input:: *)
If[UseAT,FINALPLOTAT=Plot[{FINALSAT[[2,1]]+FINALSAT[[2,2]],FINALSAT[[2,1]]-FINALSAT[[2,2]]},{t,IPat-0.5,FPat+0.5},Filling->True];
If[AutomaticAT,FullPlotAT=Show[meffATplot,FINALPLOTAT,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]]]]]]
If[Not[AutomaticAT],FullPlotAT=Show[meffATplot,FINALPLOTAT,PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT],", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=", ToString[NumberForm[m,{Digitsm,DPm}]],", \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.=",ToString[NumberForm[FinalChiSq,{3,2}]],"\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\"Manual\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)\!\(\*
StyleBox[\")\",\nFontSlant->\"Italic\",\nFontColor->RGBColor[0.5, 0, 0.5]]\)"]]]


(* ::Input:: *)
If[UseAT,Export[StringJoin[Outputpath,"/axialtensor.pdf"],FullPlotAT,"PDF"], Export[StringJoin[Outputpath,"/axialtensor.pdf"], ""]];
If[UseAT,HistogramPlotAT=Histogram[ATBoots[[;;,-1]],{BinWidth},"Probability",ImageSize->600,ChartStyle->Darker[Blue],FrameLabel->{"\!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\)/d.o.f.","Normalised Frequency"},LabelStyle->Directive[20,Black],Frame->True,PlotRange->Automatic,AxesStyle->{Thick,Black},PlotLabel->StringJoin["\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(0\)]\)\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\). Sp(",ToString[Nc],"), ",ToString[Superscript[NS,"3"],FormatType->StandardForm],"\[Times]",ToString[NT], ", \[Beta]=",ToString[NumberForm[beta,{Digitsb,DPb}]],", \!\(\*SubscriptBox[\(m\), \(o\)]\)=",ToString[NumberForm[m,{Digitsm,DPm}]]]]]


(* ::Input:: *)
If[UseAT,Export[StringJoin[Outputpath,"/axialtensorCSD.pdf"],HistogramPlotAT,"PDF"], Export[StringJoin[Outputpath,"/axialtensorCSD.pdf"], ""]];


(* ::Input:: *)
(*Square measurements, if necessary, accounting for propagation of errors and renormalise decay constants.*)
If[UseAT,massAT=Around[FINALSAT[[2,1]],FINALSAT[[2,2]]];
massAT=(GradientFlow*massAT)^2;
FINALSAT[[2,1]]=massAT["Value"];
FINALSAT[[2,2]]=massAT["Error"];]


(* ::Input:: *)
If[UseAT,Export[StringJoin[Outputpath,"/output_axialtensor.txt"],{Flatten[FINALSAT],ToString[IPat],ToString[FPat]},"Table","FieldSeparators"->"\n"], Export[StringJoin[Outputpath,"/output_axialtensor.txt"], ""]];


(* ::Section:: *)
(* Print in L AT EX form *)


(* ::Input:: *)
Export[StringJoin[Outputpath,"/",FilenameSlug,".txt"],{{"Obs","Value","Error","Chi2","Interval"},
If[UsePS,{StringJoin["\\multirow[t]{9}{*}{",ToString[NumberForm[m,{Digitsm,DPm}]],"}&$\\hat{m}_{\\text{PS}}$ &",ToString[NumberForm[Sqrt[Around[FINALSPS[[2,1]],FINALSPS[[2,2]]]]["Value"],Ceiling[Log[Accuracy[NumberForm[Sqrt[Around[FINALSPS[[2,1]],FINALSPS[[2,2]]]]["Error"],2]]]]+1]],"(",ToString[StringTake[ToString[NumberForm[Sqrt[Around[FINALSPS[[2,1]],FINALSPS[[2,2]]]]["Error"],2]],-2]],") &\\multirow{2}{*}{",ToString[NumberForm[FINALSPS[[4,1]],{3,2}]],"}& \\multirow{2}{*}{[",ToString[IPps],",",ToString[FPps],"]} \\\\"]},{"&$\\hat{m}_{\\text{PS}}]$&\\multicolumn{3}{c|}{---}\\\\"}],
If[UsePS,{"\\cline{2-3}"},{"\\cline{2-5}"}],
If[UsePS,{StringJoin["&$\\hat{f}_{\\text{PS}}$ &",ToString[NumberForm[Sqrt[Around[FINALSPS[[3,1]],FINALSPS[[3,2]]]]["Value"],Ceiling[Log[Accuracy[NumberForm[Sqrt[Around[FINALSPS[[3,1]],FINALSPS[[3,2]]]]["Error"],2]]]]+1]],"(",ToString[StringTake[ToString[NumberForm[Sqrt[Around[FINALSPS[[3,1]],FINALSPS[[3,2]]]]["Error"],2]],-2]],") && \\\\"]},{"&$\\hat{f}_{\\text{PS}}$&\\multicolumn{3}{c|}{---}\\\\"}],
{"\\cline{2-5}"},
If[UseV,{StringJoin["&$\\hat{m}_{\\text{V}}$ &",ToString[NumberForm[Sqrt[Around[FINALSV[[2,1]],FINALSV[[2,2]]]]["Value"],Ceiling[Log[Accuracy[NumberForm[Sqrt[Around[FINALSV[[2,1]],FINALSV[[2,2]]]]["Error"],2]]]]+1]],"(",ToString[StringTake[ToString[NumberForm[Sqrt[Around[FINALSV[[2,1]],FINALSV[[2,2]]]]["Error"],2]],-2]],") &\\multirow{2}{*}{",ToString[NumberForm[FINALSV[[3,1]],{3,2}]],"}& \\multirow{2}{*}{[",ToString[IPv],",",ToString[FPv],"]} \\\\"]},{"&$\\hat{m}_{\\text{V}}$&\\multicolumn{3}{c|}{---}\\\\"}],
If[UseV,{"\\cline{2-3}"},{"\\cline{2-5}"}],
If[UseV,{StringJoin["&$\\hat{f}_{\\text{V}}$ &",ToString[NumberForm[Sqrt[Around[FINALSV[[1,1]],FINALSV[[1,2]]]]["Value"],Ceiling[Log[Accuracy[NumberForm[Sqrt[Around[FINALSV[[1,1]],FINALSV[[1,2]]]]["Error"],2]]]]+1]],"(",ToString[StringTake[ToString[NumberForm[Sqrt[Around[FINALSV[[1,1]],FINALSV[[1,2]]]]["Error"],2]],-2]],") && \\\\"]},{"&$\\hat{f}_{\\text{V}}$&\\multicolumn{3}{c|}{---}\\\\"}],
{"\\cline{2-5}"},
If[UseAV,{StringJoin["&$\\hat{m}_{\\text{AV}}$ &",ToString[NumberForm[Sqrt[Around[FINALSAV[[2,1]],FINALSAV[[2,2]]]]["Value"],Ceiling[Log[Accuracy[NumberForm[Sqrt[Around[FINALSAV[[2,1]],FINALSAV[[2,2]]]]["Error"],2]]]]+1]],"(",ToString[StringTake[ToString[NumberForm[Sqrt[Around[FINALSAV[[2,1]],FINALSAV[[2,2]]]]["Error"],2]],-2]],") &\\multirow{2}{*}{",ToString[NumberForm[FINALSAV[[3,1]],{3,2}]],"}& \\multirow{2}{*}{[",ToString[IPav],",",ToString[FPav],"]} \\\\"]},{"&$\\hat{m}_{\\text{AV}}$&\\multicolumn{3}{c|}{---}\\\\"}],
If[UseAV,{"\\cline{2-3}"},{"\\cline{2-5}"}],
If[UseAV,{StringJoin["&$\\hat{f}_{\\text{AV}}$ &",ToString[NumberForm[Sqrt[Around[FINALSAV[[1,1]],FINALSAV[[1,2]]]]["Value"], Ceiling[Log[Accuracy[NumberForm[Sqrt[Around[FINALSAV[[1,1]],FINALSAV[[1,2]]]]["Error"],2]]]]+1]],"(",ToString[StringTake[ToString[NumberForm[Sqrt[Around[FINALSAV[[1,1]],FINALSAV[[1,2]]]]["Error"],2]],-2]],") && \\\\"]},{"&$\\hat{f}_{\\text{AV}}$&\\multicolumn{3}{c|}{---}\\\\"}],
{"\\cline{2-5}"},
If[UseS,{StringJoin["&$\\hat{m}_{\\text{S}}$ &",ToString[NumberForm[Sqrt[Around[FINALSS[[2,1]],FINALSS[[2,2]]]]["Value"],Ceiling[Log[Accuracy[NumberForm[Sqrt[Around[FINALSS[[2,1]],FINALSS[[2,2]]]]["Error"],2]]]]+1]],"(",ToString[StringTake[ToString[NumberForm[Sqrt[Around[FINALSS[[2,1]],FINALSS[[2,2]]]]["Error"],2]],-2]],")","&",ToString[NumberForm[FINALSS[[3,1]],{3,2}]],"& [",ToString[IPs],",",ToString[FPs],"] \\\\"]},{"&$\\hat{m}_{\\text{S}}$&\\multicolumn{3}{c|}{---}\\\\"}],
{"\\cline{2-5}"},
If[UseT,{StringJoin["&$\\hat{m}_{\\text{T}}$ &",ToString[NumberForm[Sqrt[Around[FINALST[[2,1]],FINALST[[2,2]]]]["Value"],Ceiling[Log[Accuracy[NumberForm[Sqrt[Around[FINALST[[2,1]],FINALST[[2,2]]]]["Error"],2]]]]+1]],"(",ToString[StringTake[ToString[NumberForm[Sqrt[Around[FINALST[[2,1]],FINALST[[2,2]]]]["Error"],2]],-2]],")","&",ToString[NumberForm[FINALST[[3,1]],{3,2}]],"& [",ToString[IPt],",",ToString[FPt],"] \\\\"]},{"&$\\hat{m}_{\\text{T}}$&\\multicolumn{3}{c|}{---}\\\\"}],
{"\\cline{2-5}"},
If[UseAT,{StringJoin["&$\\hat{m}_{\\text{AT}}$ &",ToString[NumberForm[Sqrt[Around[FINALSAT[[2,1]],FINALSAT[[2,2]]]]["Value"],Ceiling[Log[Accuracy[NumberForm[Sqrt[Around[FINALSAT[[2,1]],FINALSAT[[2,2]]]]["Error"],2]]]]+1]],"(",ToString[StringTake[ToString[NumberForm[Sqrt[Around[FINALSAT[[2,1]],FINALSAT[[2,2]]]]["Error"],2]],-2]],")","&",ToString[NumberForm[FINALSAT[[3,1]],{3,2}]],"& [",ToString[IPat],",",ToString[FPat],"] \\\\"]},{"&$\\hat{m}_{\\text{AT}}$&\\multicolumn{3}{c|}{---}\\\\"}],
{"\\hline\\hline"}},"Table"];


(* ::Section:: *)
(* Results and automation *)


(* ::Input:: *)
If[Not[UsePS],Print["Pseudoscalar not computed."]]
If[Not[AutomaticPS],Print["Pseudoscalar not automated."]]
If[Not[UseV],Print["Vector not computed."]]
If[Not[AutomaticV],Print["Vector not automated."]]
If[Not[UseAV],Print["Axial-vector not computed."]]
If[Not[AutomaticAV],Print["Axial-vector not automated."]]
If[Not[UseS],Print["Scalar not computed."]]
If[Not[AutomaticS],Print["Scalar not automated."]]
If[Not[UseT],Print["Tensor not computed."]]
If[Not[AutomaticT],Print["Tensor not automated."]]
If[Not[UseAT],Print["Axial-tensor not computed."]]
If[Not[AutomaticAT],Print["Axial-tensor not automated."]]
