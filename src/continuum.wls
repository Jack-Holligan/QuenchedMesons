#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Initialisation *)


(* ::Input:: *)
ClearAll;
Needs["ErrorBarPlots`"];
Needs["PlotLegends`"];
Needs["MaTeX`"];

SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1}];

Nc:=_SED_NC_;
Rep:=_SED_REP_;
Basepath:="_SED_BASEPATH_";
Outputpath:="_SED_OUTPUTPATH_";
volumes:={_SED_VOLUMES_};
Data:={"beta", "baremass", "GradientFlow", "GradientFlowErr", "Observable", "ObservableErr"};
betas:={_SED_BETAS_};
colours:={_SED_COLOURS_}; (*Blue, Orange, Darker[Green], Red, Purple *)
f[x_,L_,m2_,W_,a_]:=x*(1+L*m2)+W*a;
Observable="_SED_OBSERVABLE_";
channel="_SED_CHANNEL_";
shortchannel=<| "pseudoscalar" -> "PS", "vector" -> "V", "axialvector" -> "AV", "scalar" -> "S", "tensor" -> "T", "axialtensor" -> "AT" |>[[channel]];
SeedRandom[_SED_SEED_,Method->"MersenneTwister"];



(* ::Section:: *)
(* Computation *)

(* ::Input:: *)
Clear[Evaluate["Use"<>Observable<>channel]];
Evaluate[Symbol["Use"<>Observable<>channel]]=True;


(* ::Input:: *)
scrap={};
For[i=1,i<=Length[betas],++i,{Clear[Evaluate[StringJoin[Observable<>channel<>ToString[i]]]],If[FileExistsQ[StringJoin[Basepath,"/",ToString[Rep],"_data/",ToString[volumes[[i]]],"_",Observable,"_",channel,"_",ToString[Rep],".txt"]],Evaluate[Symbol[StringJoin[Observable<>channel<>ToString[i]]]]=Import[StringJoin[Basepath,"/",ToString[Rep],"_data/",ToString[volumes[[i]]],"_",Observable,"_",channel,"_",ToString[Rep],".txt"],"CSV","Numeric"->False],AppendTo[scrap,i]]}]


(* ::Input:: *)
Betas=Volumes=Colours={};
For[i=1,i<=Length[betas],++i,If[MemberQ[scrap,i],Continue[],{AppendTo[Betas,betas[[i]]],AppendTo[Volumes,volumes[[i]]],AppendTo[Colours,colours[[i]]]}]]
(*Remove beta values for which there is no data*)


(* ::Input:: *)
For[i=1,i<=Length[Betas],++i,{
Clear[Evaluate[StringJoin["Xdata",ToString[i]]]];
Evaluate[Symbol[StringJoin["Xdata",ToString[i]]]]=Import[StringJoin[Basepath,"/",ToString[Rep],"_data/",ToString[Volumes[[i]]],"_masses_pseudoscalar_",ToString[Rep],".txt"],"CSV","Numeric"->False]}](*Import as strings to avoid loss of trailing zeros when importing the data files.*)


(* ::Input:: *)
For[i=1,i<=Length[Betas],++i,{
Clear[Evaluate[StringJoin[Observable,channel,ToString[i]]]];
Evaluate[Symbol[StringJoin[Observable,channel,ToString[i]]]]=Import[StringJoin[Basepath,"/",ToString[Rep],"_data/",ToString[Volumes[[i]]],"_",Observable,"_",channel,"_",ToString[Rep],".txt"],"CSV","Numeric"->False]}]


(* ::Input:: *)
For[i=1,i<=Length[Betas],++i,Clear[Evaluate["massPS"<>ToString[i]]];tmp={};
For[j=1,j<=Dimensions[Symbol[StringJoin[Observable,channel,ToString[i]]]][[1]],++j,
mass=Symbol[StringJoin[Observable,channel,ToString[i]]][[j,2]];
place=Position[Symbol[Evaluate["Xdata"<>ToString[i]]],mass][[1,1]];
AppendTo[tmp,Evaluate[Symbol["Xdata"<>ToString[i]]][[place]]]];
Evaluate[Symbol["massPS"<>ToString[i]]]=tmp]
(*Include only those bare masses in massPS that are also in the observable's array*)


(* ::Input:: *)
massPSboots={};For[i=1,i<=Length[Betas],++i,For[j=1,j<=Length[Evaluate[Symbol["massPS"<>ToString[i]]][[;;,2]]],++j,
{row=Flatten[Import[StringJoin[Basepath,"/",ToString[Rep],"_data/",ToString[Volumes[[i]]],"_m",ToString[Rep],ToString[Evaluate[Symbol["massPS"<>ToString[i]][[j,2]]]],"_pseudoscalar_masses_boots.csv"]]],gradient=ToExpression[Evaluate[Symbol["massPS"<>ToString[i]]][[j,3]]],AppendTo[massPSboots,gradient^2*row]}]];
(*Import bootstrap data for pseudoscalar mass multiplying by the relevant gradient flow scale*)


(* ::Input:: *)
Clear[Evaluate[StringJoin[Observable,channel,"boots"]]];
tmpboots={};For[i=1,i<=Length[Betas],++i,For[j=1,j<=Length[Evaluate[Symbol[StringJoin[Observable,channel,ToString[i]]]][[;;,2]]],++j,{row=Flatten[Import[StringJoin[Basepath,"/",ToString[Rep],"_data/",ToString[Volumes[[i]]],"_m",ToString[Rep],ToString[Evaluate[Symbol["massPS"<>ToString[i]][[j,2]]]],"_",channel,"_",Observable,"_boots.csv"]]],gradient=ToExpression[Evaluate[Symbol[Observable<>channel<>ToString[i]]][[j,3]]],
AppendTo[tmpboots,gradient^2*row]}]];
Evaluate[Symbol[Observable<>channel<>"boots"]]=tmpboots;(*Import bootstrap data for observable multiplying by the relevant gradient flow scale*)


(* ::Input:: *)
For[i=1,i<=Length[Betas],++i,{holder=Symbol[Evaluate[Observable<>channel<>ToString[i]]],
For[j=1,j<=Dimensions[holder][[1]],++j,
holder[[j]]=Flatten@Map[Interpreter["Number"],StringSplit[holder[[j]]]]],Clear[Evaluate[Observable<>channel<>ToString[i]]],Evaluate[Symbol[Observable<>channel<>ToString[i]]]=holder}]
For[i=1,i<=Length[Betas],++i,{holder=Symbol[Evaluate["massPS"<>ToString[i]]],
For[j=1,j<=Dimensions[holder][[1]],++j,
holder[[j]]=Flatten@Map[Interpreter["Number"],StringSplit[holder[[j]]]]],Clear[Evaluate["massPS"<>ToString[i]]],Evaluate[Symbol["massPS"<>ToString[i]]]=holder}]


(* ::Input:: *)
For[i=1,i<=Length[Betas],++i,{For[j=1,j<=Dimensions[Evaluate[Symbol[Observable<>channel<>ToString[i]]]][[1]],++j,{
holder1=Symbol["massPS"<>ToString[i]];
holder2=Symbol[Observable<>channel<>ToString[i]];
gradient=Around[Symbol["massPS"<>ToString[i]][[1,3]],Symbol["massPS"<>ToString[i]][[1,4]]];
tmp1=Around[holder1[[j,5]],holder1[[j,6]]]*gradient^2;
tmp2=Around[holder2[[j,5]],holder2[[j,6]]]*gradient^2;
holder1[[j,5]]=tmp1["Value"];
holder1[[j,6]]=tmp1["Error"];
holder2[[j,5]]=tmp2["Value"];
holder2[[j,6]]=tmp2["Error"];
Clear[Evaluate["massPS"<>ToString[i]]];Evaluate[Symbol["massPS"<>ToString[i]]]=holder1;
Clear[Evaluate[Observable<>channel<>ToString[i]]];Evaluate[Symbol[Observable<>channel<>ToString[i]]]=holder2;}]}](*Multiply all observables by the relevant gradient flow scale.*)


(* ::Input:: *)
ChannelSymbol=If[ToString[Rep]=="S","\\mathcal{"<>shortchannel<>"}","\\mathrm{"<>If[ToString[Rep]=="AS",ToLowerCase[shortchannel],shortchannel]<>"}"];
If[StringMatchQ[Observable,"decayconsts"],ylabel=MaTeX["\\hat{f}^2_{"<>ChannelSymbol<>"}",FontSize->24],If[StringMatchQ[Observable,"masses"],ylabel=MaTeX["\\hat{m}^2_{"<>ChannelSymbol<>"}",FontSize->24]]];
xvalues={};yvalues={};
For[i=1,i<=Length[Betas],++i,{AppendTo[xvalues,Evaluate[Symbol["massPS"<>ToString[i]]][[;;,5]]+Evaluate[Symbol["massPS"<>ToString[i]]][[;;,6]]],AppendTo[xvalues,Evaluate[Symbol["massPS"<>ToString[i]]][[;;,5]]-Evaluate[Symbol["massPS"<>ToString[i]]][[;;,6]]]}]
For[i=1,i<=Length[Betas],++i,{AppendTo[yvalues,Evaluate[Symbol[Observable<>channel<>ToString[i]]][[;;,5]]+Evaluate[Symbol[Observable<>channel<>ToString[i]]][[;;,6]]],AppendTo[yvalues,Evaluate[Symbol[Observable<>channel<>ToString[i]]][[;;,5]]-Evaluate[Symbol[Observable<>channel<>ToString[i]]][[;;,6]]]}]


(* ::Input:: *)
maxmassPS=_SED_MAXMASSPS_;
minmassPS=_SED_MINMASSPS_;


(* ::Input:: *)
Clear[Evaluate["dof"<>Observable<>channel]];doft=0;
For[i=1,i<=Length[Betas],++i,doft+=Length[Select[Select[Evaluate[Symbol["massPS"<>ToString[i]]][[;;,5]],#<=maxmassPS&],#>=minmassPS&]]]
Evaluate[Symbol["dof"<>Observable<>channel]]=doft-3;


(* ::Input:: *)
tmp=NMinimize[{Sum[Sum[Boole[Evaluate[Symbol["massPS"<>ToString[i]]][[k,5]]<=maxmassPS&&Evaluate[Symbol["massPS"<>ToString[i]]][[k,5]]>=minmassPS]*(f[x,L,Evaluate[Symbol["massPS"<>ToString[i]]][[k,5]],W,1/Evaluate[Symbol["massPS"<>ToString[i]]][[k,3]]]-Evaluate[Symbol[Observable<>channel<>ToString[i]][[k,5]]])^2/(Evaluate[Symbol[Observable<>channel<>ToString[i]][[k,6]]])^2,{k,1,Dimensions[Evaluate[Symbol["massPS"<>ToString[i]]]][[1]]}],{i,1,Length[Betas]}],x>=0},{x,L,W},Method->"DifferentialEvolution"];
Clear[Evaluate[Observable<>channel<>"final"]];
Evaluate[Symbol[Observable<>channel<>"final"]]={x/.tmp[[2,1]], L/.tmp[[2,2]], W/.tmp[[2,3]], tmp[[1]]/Evaluate[Symbol["dof"<>Observable<>channel]]}


(* ::Input:: *)
For[i=1,i<=Length[Betas],++i,
Clear[Evaluate["plot"<>Observable<>channel<>ToString[i]]];
Evaluate[Symbol["plot"<>Observable<>channel<>ToString[i]]]=ErrorListPlot[Table[{{Evaluate[Symbol["massPS"<>ToString[i]]][[t,5]],Evaluate[Symbol[StringJoin[Observable,channel,ToString[i]]]][[t,5]]},ErrorBar[Evaluate[Symbol["massPS"<>ToString[i]]][[t,6]],Evaluate[Symbol[StringJoin[Observable,channel,ToString[i]]]][[t,6]]]},{t,Dimensions[Evaluate[Symbol["massPS"<>ToString[i]]]][[1]]}],Frame->True,FrameLabel->{MaTeX["\\hat{m}^2_{\\rm PS}",FontSize->24],ylabel},LabelStyle->Directive[18,Black,FontFamily->Times],RotateLabel->False,PlotStyle->Echo[Colours[[i]]] ,PlotLegends->{MaTeX["\\beta="<>ToString[NumberForm[Betas[[i]],{3,2}]],FontSize->24]},PlotRange->{{-0.1,1.1Max[xvalues]},{0,1.1Max[yvalues]}}]];


(* ::Input:: *)
showplots={};
For[i=1,i<=Length[Betas],++i,AppendTo[showplots,Symbol["plot"<>Observable<>channel<>ToString[i]]]];
Show[showplots]


(* ::Input:: *)
(*Export[Outputpath<>"/"<>ToString[Rep]<>"/"<>channel<>"_"<>Observable<>"_"<>ToString[Rep]<>"_Sp"<>ToString[Nc]<>".pdf",Show[showplots],"PDF"];*)


(* ::Input:: *)
Nboots=Dimensions[massPSboots][[2]];


(* ::Input:: *)
massPSave=Table[Mean[massPSboots[[i,;;]]],{i,1,Length[massPSboots[[;;,1]]]}];
Clear[Evaluate[Observable<>channel<>"ave"]];Evaluate[Symbol[Observable<>channel<>"ave"]]=Table[Mean[Evaluate[Symbol[Observable<>channel<>"boots"]][[i,;;]]],{i,1,Length[massPSboots[[;;,1]]]}];
(*Compute the average of the PSmass and observable's bootstrap samples.*)


(* ::Input:: *)
Clear[Evaluate[StringJoin[Observable,channel,"cov"]]];
Evaluate[Symbol[StringJoin[Observable,channel,"cov"]]]=
Table[Table[Mean[Table[(Evaluate[Symbol[Observable<>channel<>"boots"]][[j,i]]-Evaluate[Symbol[Observable<>channel<>"ave"]][[j]])(Evaluate[Symbol[Observable<>channel<>"boots"]][[k,i]]-Evaluate[Symbol[StringJoin[Observable,channel,"ave"]]][[k]]),{i,1,Nboots}]],{j,1,Length[Evaluate[Symbol[Observable<>channel<>"ave"]]]}],{k,1,Length[Evaluate[Symbol[Observable<>channel<>"ave"]]]}];(*Compute the covariance matrix for the observable's data.*)


(* ::Input:: *)
Clear[Evaluate["InvCov"<>Observable<>channel]];
Evaluate[Symbol["InvCov"<>Observable<>channel]]=Inverse[Evaluate[Symbol[Observable<>channel<>"cov"]]];


(* ::Input:: *)
BootstrappingXdata={};BootstrappingYdata={};
For[i=1,i<=Length[Betas],++i,{BootstrappingXdata=Join[BootstrappingXdata,Evaluate[Symbol["massPS"<>ToString[i]]]],BootstrappingYdata=Join[BootstrappingYdata,Evaluate[Symbol[Observable<>channel<>ToString[i]]]]}]


(* ::Input:: *)
Clear[x,L,W]
ProgressIndicator[Dynamic[n],{1,Nboots}]


(* ::Input:: *)
Clear[Evaluate[channel<>Observable<>"results"]]
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{Evaluate[Symbol[channel<>Observable<>"results"]]=Table[n=k;tmp=NMinimize[{Sum[Sum[(f[x,L,massPSboots[[i,k]],W,1/BootstrappingXdata[[i,3]]]-Evaluate[Symbol[Observable<>channel<>"boots"]][[i,k]])Evaluate[Symbol["InvCov"<>Observable<>channel]][[i,j]](f[x,L,massPSboots[[j,k]],W,1/BootstrappingXdata[[j,3]]]-Evaluate[Symbol[Observable<>channel<>"boots"]][[j,k]]),{i,1,Length[BootstrappingXdata[[;;,1]]]}],{j,1,Length[BootstrappingXdata[[;;,1]]]}],0.8*Evaluate[Symbol[Observable<>channel<>"final"]][[1]]<=x<=1.2*Evaluate[Symbol[Observable<>channel<>"final"]][[1]]},{x,L,W},Method->"NelderMead"];{x/.tmp[[2,1]],L/.tmp[[2,2]],W/.tmp[[2,3]],tmp[[1]]/Evaluate[Symbol["dof"<>Observable<>channel]]},{k,1,Nboots}];}];


(* ::Input:: *)
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{resultsholder={};}];


(* ::Input:: *)
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{data=Abs[Evaluate[Symbol[channel<>Observable<>"results"]][[;;,1]]];
gaussian=DistributionFitTest[data,Automatic,"HypothesisTestData"]["FittedDistribution"];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[gaussian,x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[gaussian[[1]]],",",ToString[gaussian[[2]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{channel<>" "<>Observable<>", x","PDF"},LabelStyle->Directive[15,Black],Frame->True]};]
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{AppendTo[resultsholder,{gaussian[[1]],gaussian[[2]]}];}];


(* ::Input:: *)
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{data=Abs[Evaluate[Symbol[channel<>Observable<>"results"]][[;;,2]]];}];
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{gaussian=DistributionFitTest[data,Automatic,"HypothesisTestData"]["FittedDistribution"];}];
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{Show[Histogram[data,Automatic,"PDF"],Plot[PDF[gaussian,x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[gaussian[[1]]],",",ToString[gaussian[[2]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{channel<>" "<>Observable<>", L","PDF"},LabelStyle->Directive[15,Black],Frame->True]}];
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{AppendTo[resultsholder,{gaussian[[1]],gaussian[[2]]}];}];


(* ::Input:: *)
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{data=Abs[Evaluate[Symbol[channel<>Observable<>"results"]][[;;,3]]];}];
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{gaussian=DistributionFitTest[data,Automatic,"HypothesisTestData"]["FittedDistribution"];}];
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{Show[Histogram[data,Automatic,"PDF"],Plot[PDF[gaussian,x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[gaussian[[1]]],",",ToString[gaussian[[2]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{channel<>" "<>Observable<>", W","PDF"},LabelStyle->Directive[15,Black],Frame->True]}];
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{AppendTo[resultsholder,{gaussian[[1]],gaussian[[2]]}];}];


(* ::Input:: *)
Clear[Evaluate["FINALS"<>channel<>Observable]];
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{Evaluate[Symbol["FINALS"<>channel<>Observable]]=ArrayReshape[Riffle[Evaluate[Symbol[Observable<>channel<>"final"]],resultsholder[[;;,2]]],{4,2}]}];


(* ::Input:: *)
Clear[Evaluate["dash"<>channel<>Observable],Evaluate["limit"<>channel<>"1"],Evaluate["limit"<>channel<>"2"]]
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{Evaluate[Symbol["dash"<>channel<>Observable]]=ErrorListPlot[{{{0,Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]]},ErrorBar[Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]]}},PlotMarkers->{"\[FilledSquare]",15},PlotRange->{{-0.1,1.1Max[xvalues]},{0,1.1Max[yvalues,Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]]+Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]]}},AxesOrigin->{-0.1,0},Frame->True,FrameLabel->{"\!\(\*SubsuperscriptBox[OverscriptBox[\(m\), \(^\)], \(PS\), \(2\)]\)",ylabel},LabelStyle->Directive[15,Black],PlotLegends->{MaTeX["\\text{Extrapolation}",FontSize->24]}];
AppendTo[showplots,Evaluate[Symbol["dash"<>channel<>Observable]]];
line=Plot[{f[Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]],Evaluate[Symbol["FINALS"<>channel<>Observable]][[2,1]],x,Evaluate[Symbol["FINALS"<>channel<>Observable]][[3,1]],0]+Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]],f[Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]],Evaluate[Symbol["FINALS"<>channel<>Observable]][[2,1]],x,Evaluate[Symbol["FINALS"<>channel<>Observable]][[3,1]],0]-Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]},{x,0,1.1Max[xvalues]},PlotStyle->{{Dashed,Black,Thickness[0.001]},{Dashed,Black,Thickness[0.001]}},Filling->True,AxesOrigin->{-0.1,0},Frame->True,FrameLabel->{"\!\(\*SubsuperscriptBox[OverscriptBox[\(m\), \(^\)], \(PS\), \(2\)]\)",ylabel},LabelStyle->Directive[22,Black]];
AppendTo[showplots,line];
Evaluate[Symbol["limit"<>channel<>"1"]]=Rectangle[{maxmassPS,0},{1.1Max[xvalues],1.1Max[yvalues,Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]]+Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]]}];
Evaluate[Symbol["limit"<>channel<>"2"]]=Rectangle[{If[minmassPS!=0,Min[0.1xvalues],0],0},{minmassPS,1.1Max[yvalues,Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]]+Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]]}];}];


(* ::Input:: *)
Clear[Evaluate["dash"<>channel<>Observable],Evaluate["limit"<>channel<>"1"],Evaluate["limit"<>channel<>"2"]]
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{Evaluate[Symbol["dash"<>channel<>Observable]]=ErrorListPlot[{{{0,Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]]},ErrorBar[Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]]}},PlotMarkers->{"\[FilledUpTriangle]",15},PlotStyle->Black,PlotRange->{{-0.1,1.1Max[xvalues]},{0,1.1Max[yvalues,Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]]+Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]]}},AxesOrigin->{-0.1,0},Frame->True,FrameLabel->{"\!\(\*SubsuperscriptBox[OverscriptBox[\(m\), \(^\)], \(PS\), \(2\)]\)",ylabel},LabelStyle->Directive[15,Black],PlotLegends->{MaTeX["\\text{Extrapolation (higher }\\beta)",FontSize->24]}];
AppendTo[showplots,Evaluate[Symbol["dash"<>channel<>Observable]]];
line=Plot[{f[Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]],Evaluate[Symbol["FINALS"<>channel<>Observable]][[2,1]],x,Evaluate[Symbol["FINALS"<>channel<>Observable]][[3,1]],0]+Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]],f[Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]],Evaluate[Symbol["FINALS"<>channel<>Observable]][[2,1]],x,Evaluate[Symbol["FINALS"<>channel<>Observable]][[3,1]],0]-Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]},{x,0,1.1Max[xvalues]},PlotStyle->{{Dashed,Gray,Opacity[0.0001]},{Dashed,Gray,Opacity[0.0001]}},Filling->True,AxesOrigin->{-0.1,0},Frame->True,FrameLabel->{"\!\(\*SubsuperscriptBox[OverscriptBox[\(m\), \(^\)], \(PS\), \(2\)]\)",ylabel},LabelStyle->Directive[22,Black]];
AppendTo[showplots,line];
Evaluate[Symbol["limit"<>channel<>"1"]]=Rectangle[{maxmassPS,0},{1.1Max[xvalues],1.1Max[yvalues,Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]]+Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]]}];
Evaluate[Symbol["limit"<>channel<>"2"]]=Rectangle[{If[minmassPS!=0,Min[0.1xvalues],0],0},{minmassPS,1.1Max[yvalues,Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,1]]+Evaluate[Symbol["FINALS"<>channel<>Observable]][[1,2]]]}];}];


(* ::Input:: *)
FitBand=Show[Evaluate[Symbol["dash"<>channel<>Observable]],line]


(* ::Input:: *)
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,Show[showplots],PlotLabel->MaTeX["\\chi^2/\\text{d.o.f.}="<>ToString[NumberForm[Evaluate[Symbol[Observable<>channel<>"final"]][[4]],{3,2}]],FontSize->24],Epilog->{Directive[Red,Opacity[0.1]],Evaluate[Symbol["limit"<>channel<>"1"]],Evaluate[Symbol["limit"<>channel<>"2"]]}]]


(* ::Input:: *)
FinalPlot=If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,Show[showplots,PlotLabel->MaTeX["\\chi^2/\\text{d.o.f.}="<>ToString[NumberForm[Evaluate[Symbol[Observable<>channel<>"final"]][[4]],{3,2}]],FontSize->24],Epilog->{Directive[Red,Opacity[0.1]],Evaluate[Symbol["limit"<>channel<>"1"]],Evaluate[Symbol["limit"<>channel<>"2"]]}]]


(* ::Input:: *)
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{Export[Outputpath<>"/"<>ToString[Rep]<>"/"<>channel<>"_"<>Observable<>"_"<>ToString[Rep]<>"_Sp"<>ToString[Nc]<>".pdf",FinalPlot,"PDF",AllowRasterization->False,ImageResolution->1200];}];


(* ::Input:: *)
If[Evaluate[Symbol["Use"<>Observable<>channel]]==True,{Export[Outputpath<>"/"<>ToString[Rep]<>"/"<>channel<>"_"<>Observable<>"_"<>ToString[Rep]<>"_Sp"<>ToString[Nc]<>".dat",Evaluate[Symbol["FINALS"<>channel<>Observable]],"Table","FieldSeparators"->" "];},{Export[Basepath<>"/"<>ToString[Rep]<>"/"<>channel<>"_"<>Observable<>"_"<>ToString[Rep]<>"_Sp"<>ToString[Nc]".dat",Array["--"&,{4,2}],"Table","FieldSeparators"->" "];}];
