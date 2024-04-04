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

OutputFilenameBase=Outputpath<>"/"<>ToString[Rep]<>"/"<>channel<>"_"<>Observable<>"_"<>ToString[Rep]<>"_Sp"<>ToString[Nc];


(* ::Section:: *)
(* Computation *)

(* ::Input:: *)
UseObservable=_SED_USE_;


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
ObservableBoots={};For[i=1,i<=Length[Betas],++i,For[j=1,j<=Length[Evaluate[Symbol[StringJoin[Observable,channel,ToString[i]]]][[;;,2]]],++j,{row=Flatten[Import[StringJoin[Basepath,"/",ToString[Rep],"_data/",ToString[Volumes[[i]]],"_m",ToString[Rep],ToString[Evaluate[Symbol["massPS"<>ToString[i]][[j,2]]]],"_",channel,"_",Observable,"_boots.csv"]]],gradient=ToExpression[Evaluate[Symbol[Observable<>channel<>ToString[i]]][[j,3]]],
AppendTo[ObservableBoots,gradient^2*row]}]];
(*Import bootstrap data for observable multiplying by the relevant gradient flow scale*)


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
ChannelSymbol[shortchannel_]:=If[ToString[Rep]=="S","\\mathcal{"<>shortchannel<>"}","\\mathrm{"<>If[ToString[Rep]=="AS",ToLowerCase[shortchannel],shortchannel]<>"}"];
If[StringMatchQ[Observable,"decayconsts"],ylabel=MaTeX["\\hat{f}^2_{"<>ChannelSymbol[shortchannel]<>"}",FontSize->20],If[StringMatchQ[Observable,"masses"],ylabel=MaTeX["\\hat{m}^2_{"<>ChannelSymbol[shortchannel]<>"}",FontSize->20]]];
xvalues={};yvalues={};
For[i=1,i<=Length[Betas],++i,{AppendTo[xvalues,Evaluate[Symbol["massPS"<>ToString[i]]][[;;,5]]+Evaluate[Symbol["massPS"<>ToString[i]]][[;;,6]]],AppendTo[xvalues,Evaluate[Symbol["massPS"<>ToString[i]]][[;;,5]]-Evaluate[Symbol["massPS"<>ToString[i]]][[;;,6]]]}]
For[i=1,i<=Length[Betas],++i,{AppendTo[yvalues,Evaluate[Symbol[Observable<>channel<>ToString[i]]][[;;,5]]+Evaluate[Symbol[Observable<>channel<>ToString[i]]][[;;,6]]],AppendTo[yvalues,Evaluate[Symbol[Observable<>channel<>ToString[i]]][[;;,5]]-Evaluate[Symbol[Observable<>channel<>ToString[i]]][[;;,6]]]}]


(* ::Input:: *)
maxmassPS=_SED_MAXMASSPS_;
minmassPS=_SED_MINMASSPS_;


(* ::Input:: *)
doft=0;
For[i=1,i<=Length[Betas],++i,doft+=Length[Select[Select[Evaluate[Symbol["massPS"<>ToString[i]]][[;;,5]],#<=maxmassPS&],#>=minmassPS&]]]
ObservableDof=doft-3;


(* ::Input:: *)
tmp=NMinimize[{Sum[Sum[Boole[Evaluate[Symbol["massPS"<>ToString[i]]][[k,5]]<=maxmassPS&&Evaluate[Symbol["massPS"<>ToString[i]]][[k,5]]>=minmassPS]*(f[x,L,Evaluate[Symbol["massPS"<>ToString[i]]][[k,5]],W,1/Evaluate[Symbol["massPS"<>ToString[i]]][[k,3]]]-Evaluate[Symbol[Observable<>channel<>ToString[i]][[k,5]]])^2/(Evaluate[Symbol[Observable<>channel<>ToString[i]][[k,6]]])^2,{k,1,Dimensions[Evaluate[Symbol["massPS"<>ToString[i]]]][[1]]}],{i,1,Length[Betas]}],x>=0},{x,L,W},Method->"DifferentialEvolution"];
ObservableFinal={x/.tmp[[2,1]], L/.tmp[[2,2]], W/.tmp[[2,3]], tmp[[1]]/ObservableDof}


(* ::Input:: *)
For[i=1,i<=Length[Betas],++i,
Evaluate[Symbol["plot"<>Observable<>channel<>ToString[i]]]=ErrorListPlot[Table[{{Evaluate[Symbol["massPS"<>ToString[i]]][[t,5]],Evaluate[Symbol[StringJoin[Observable,channel,ToString[i]]]][[t,5]]},ErrorBar[Evaluate[Symbol["massPS"<>ToString[i]]][[t,6]],Evaluate[Symbol[StringJoin[Observable,channel,ToString[i]]]][[t,6]]]},{t,Dimensions[Evaluate[Symbol["massPS"<>ToString[i]]]][[1]]}],Axes->False,Frame->True,FrameLabel->{MaTeX["\\hat{m}^2_{"<>ChannelSymbol["PS"]<>"}",FontSize->20],ylabel},LabelStyle->Directive[18,Black,FontFamily->Times],RotateLabel->False,PlotStyle->Colours[[i]] ,PlotLegends->{MaTeX["\\beta="<>ToString[NumberForm[Betas[[i]],{3,2}]],FontSize->20]},PlotRange->{{-0.1,1.1Max[xvalues]},{0,1.1Max[yvalues]}}]];


(* ::Input:: *)
showplots={};
For[i=1,i<=Length[Betas],++i,AppendTo[showplots,Symbol["plot"<>Observable<>channel<>ToString[i]]]];


(* ::Input:: *)
(*Export[OutputFilenameBase<>".pdf",Show[showplots],"PDF"];*)


(* ::Input:: *)
Nboots=Dimensions[massPSboots][[2]];


(* ::Input:: *)
ave=Table[Mean[ObservableBoots[[i,;;]]],{i,1,Length[massPSboots[[;;,1]]]}];
(*Compute the average of the PSmass and observable's bootstrap samples.*)


(* ::Input:: *)
ObservableCovariance=Table[Table[Mean[Table[(ObservableBoots[[j,i]]-ave[[j]])(ObservableBoots[[k,i]]-ave[[k]]),{i,1,Nboots}]],{j,1,Length[ave]}],{k,1,Length[ave]}];(*Compute the covariance matrix for the observable's data.*)


(* ::Input:: *)
InvCov=Inverse[ObservableCovariance];


(* ::Input:: *)
BootstrappingXdata={};BootstrappingYdata={};
For[i=1,i<=Length[Betas],++i,{BootstrappingXdata=Join[BootstrappingXdata,Evaluate[Symbol["massPS"<>ToString[i]]]],BootstrappingYdata=Join[BootstrappingYdata,Evaluate[Symbol[Observable<>channel<>ToString[i]]]]}]


(* ::Input:: *)
Clear[x,L,W]
ProgressIndicator[Dynamic[n],{1,Nboots}]


(* ::Input:: *)
If[UseObservable,{ObservableResults=Table[n=k;tmp=NMinimize[{Sum[Sum[(f[x,L,massPSboots[[i,k]],W,1/BootstrappingXdata[[i,3]]]-ObservableBoots[[i,k]])InvCov[[i,j]](f[x,L,massPSboots[[j,k]],W,1/BootstrappingXdata[[j,3]]]-ObservableBoots[[j,k]]),{i,1,Length[BootstrappingXdata[[;;,1]]]}],{j,1,Length[BootstrappingXdata[[;;,1]]]}],0.8*ObservableFinal[[1]]<=x<=1.2*ObservableFinal[[1]]},{x,L,W},Method->"NelderMead"];{x/.tmp[[2,1]],L/.tmp[[2,2]],W/.tmp[[2,3]],tmp[[1]]/ObservableDof},{k,1,Nboots}];}];


(* ::Input:: *)
If[UseObservable,{resultsholder={};}];


(* ::Input:: *)
If[UseObservable,{data=Abs[ObservableResults[[;;,1]]];
gaussian=DistributionFitTest[data,Automatic,"HypothesisTestData"]["FittedDistribution"];
Show[Histogram[data,Automatic,"PDF"],Plot[PDF[gaussian,x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[gaussian[[1]]],",",ToString[gaussian[[2]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{channel<>" "<>Observable<>", x","PDF"},LabelStyle->Directive[15,Black],Frame->True]};]
If[UseObservable,{AppendTo[resultsholder,{gaussian[[1]],gaussian[[2]]}];}];


(* ::Input:: *)
If[UseObservable,{data=Abs[ObservableResults[[;;,2]]];}];
If[UseObservable,{gaussian=DistributionFitTest[data,Automatic,"HypothesisTestData"]["FittedDistribution"];}];
If[UseObservable,{Show[Histogram[data,Automatic,"PDF"],Plot[PDF[gaussian,x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[gaussian[[1]]],",",ToString[gaussian[[2]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{channel<>" "<>Observable<>", L","PDF"},LabelStyle->Directive[15,Black],Frame->True]}];
If[UseObservable,{AppendTo[resultsholder,{gaussian[[1]],gaussian[[2]]}];}];


(* ::Input:: *)
If[UseObservable,{data=Abs[ObservableResults[[;;,3]]];}];
If[UseObservable,{gaussian=DistributionFitTest[data,Automatic,"HypothesisTestData"]["FittedDistribution"];}];
If[UseObservable,{Show[Histogram[data,Automatic,"PDF"],Plot[PDF[gaussian,x],{x,Min[data],Max[data]},PlotLegends->StringJoin["(\[Mu],\[Sigma])=(",ToString[gaussian[[1]]],",",ToString[gaussian[[2]]],")"],PlotStyle->Thick,PlotRange->Full],FrameLabel->{channel<>" "<>Observable<>", W","PDF"},LabelStyle->Directive[15,Black],Frame->True]}];
If[UseObservable,{AppendTo[resultsholder,{gaussian[[1]],gaussian[[2]]}];}];


(* ::Input:: *)
If[UseObservable,{CombinedResultsFinal=ArrayReshape[Riffle[ObservableFinal,resultsholder[[;;,2]]],{4,2}]}];


(* ::Input:: *)
If[UseObservable,{DashAllBeta=ErrorListPlot[{{{0,CombinedResultsFinal[[1,1]]},ErrorBar[CombinedResultsFinal[[1,2]]]}},PlotMarkers->{"\[FilledSquare]",15},PlotRange->{{-0.1,1.1Max[xvalues]},{0,1.1Max[yvalues,CombinedResultsFinal[[1,1]]+CombinedResultsFinal[[1,2]]]}},AxesOrigin->{-0.1,0},Frame->True,FrameLabel->{"\!\(\*SubsuperscriptBox[OverscriptBox[\(m\), \(^\)], \(PS\), \(2\)]\)",ylabel},LabelStyle->Directive[15,Black],PlotLegends->{MaTeX["\\text{Extrapolation}",FontSize->20]}];
AppendTo[showplots,DashAllBeta];
line=Plot[{f[CombinedResultsFinal[[1,1]],CombinedResultsFinal[[2,1]],x,CombinedResultsFinal[[3,1]],0]+CombinedResultsFinal[[1,2]],f[CombinedResultsFinal[[1,1]],CombinedResultsFinal[[2,1]],x,CombinedResultsFinal[[3,1]],0]-CombinedResultsFinal[[1,2]]},{x,0,1.1Max[xvalues]},PlotStyle->{{Dashed,Black,Thickness[0.001]},{Dashed,Black,Thickness[0.001]}},Filling->True,AxesOrigin->{-0.1,0},Frame->True,FrameLabel->{"\!\(\*SubsuperscriptBox[OverscriptBox[\(m\), \(^\)], \(PS\), \(2\)]\)",ylabel},LabelStyle->Directive[22,Black]];
AppendTo[showplots,line];
limit1=Rectangle[{maxmassPS,0},{1.1Max[xvalues],1.1Max[yvalues,CombinedResultsFinal[[1,1]]+CombinedResultsFinal[[1,2]]]}];
limit2=Rectangle[{If[minmassPS!=0,Min[0.1xvalues],0],0},{minmassPS,1.1Max[yvalues,CombinedResultsFinal[[1,1]]+CombinedResultsFinal[[1,2]]]}];}];


(* ::Input:: *)
If[And[Not[ContainsAll[scrap,{1}]],Length[Betas]>1,UseObservable],
doft=0;
For[i=2,i<=Length[Betas],++i,doft+=Length[Select[Select[Evaluate[Symbol["massPS"<>ToString[i]]][[;;,5]],#<=maxmassPS&],#>=minmassPS&]]];
ObservableDofHighBeta=doft-3;
SecondBetaIndex=1+Length[Evaluate[Symbol[Observable<>channel<>"1"]]];
tmp=NMinimize[{Sum[Sum[Boole[Evaluate[Symbol["massPS"<>ToString[i]]][[k,5]]<=maxmassPS&&Evaluate[Symbol["massPS"<>ToString[i]]][[k,5]]>=minmassPS]*(f[x,L,Evaluate[Symbol["massPS"<>ToString[i]]][[k,5]],W,1/Evaluate[Symbol["massPS"<>ToString[i]]][[k,3]]]-Evaluate[Symbol[Observable<>channel<>ToString[i]][[k,5]]])^2/(Evaluate[Symbol[Observable<>channel<>ToString[i]][[k,6]]])^2,{k,1,Dimensions[Evaluate[Symbol["massPS"<>ToString[i]]]][[1]]}],{i,2,Length[Betas]}],x>=0},{x,L,W},Method->"DifferentialEvolution"];
FinalHighBeta={x/.tmp[[2,1]], L/.tmp[[2,2]], W/.tmp[[2,3]], tmp[[1]]/ObservableDofHighBeta};
massPSbootsHighBeta=Table[massPSboots[[i]], {i, SecondBetaIndex, Length[massPSboots]}];
ObservableBootsHighBeta=Table[ObservableBoots[[i]], {i,SecondBetaIndex,Length[ObservableBoots]}];
covHighBeta=Table[Table[Mean[Table[(ObservableBoots[[j,i]]-ave[[j]])(ObservableBoots[[k,i]]-ave[[k]]),{i,1,Nboots}]],{j,SecondBetaIndex,Length[ave]}],{k,SecondBetaIndex,Length[ave]}];
InvCovHighBeta=Inverse[covHighBeta];
BootstrappingXdataHighBeta={};BootstrappingYdataHighBeta={};
For[i=2,i<=Length[Betas],++i,{BootstrappingXdataHighBeta=Join[BootstrappingXdataHighBeta,Evaluate[Symbol["massPS"<>ToString[i]]]],BootstrappingYdataHighBeta=Join[BootstrappingYdataHighBeta,Evaluate[Symbol[Observable<>channel<>ToString[i]]]]}];
ResultsHighBeta=Table[n=k;tmp=NMinimize[{Sum[Sum[(f[x,L,massPSbootsHighBeta[[i,k]],W,1/BootstrappingXdataHighBeta[[i,3]]]-ObservableBootsHighBeta[[i,k]])InvCovHighBeta[[i,j]](f[x,L,massPSbootsHighBeta[[j,k]],W,1/BootstrappingXdataHighBeta[[j,3]]]-ObservableBootsHighBeta[[j,k]]),{i,1,Length[BootstrappingXdataHighBeta[[;;,1]]]}],{j,1,Length[BootstrappingXdataHighBeta[[;;,1]]]}],0.8*FinalHighBeta[[1]]<=x<=1.2*FinalHighBeta[[1]]},{x,L,W},Method->"NelderMead"];{x/.tmp[[2,1]],L/.tmp[[2,2]],W/.tmp[[2,3]],tmp[[1]]/ObservableDofHighBeta},{k,1,Nboots}];
resultsholder=Table[gaussian=DistributionFitTest[Abs[ResultsHighBeta[[;;,i]]], Automatic,"HypothesisTestData"]["FittedDistribution"];{gaussian[[1]], gaussian[[2]]}, {i, 1, 3}];
CombinedResultsHighBeta=ArrayReshape[Riffle[FinalHighBeta,resultsholder[[;;,2]]],{4,2}];
DashHighBeta=ErrorListPlot[{{{0,CombinedResultsHighBeta[[1,1]]},ErrorBar[CombinedResultsHighBeta[[1,2]]]}},PlotMarkers->{"\[FilledUpTriangle]",15},PlotStyle->Black,PlotRange->{{-0.1,1.1Max[xvalues]},{0,1.1Max[yvalues,CombinedResultsHighBeta[[1,1]]+CombinedResultsHighBeta[[1,2]]]}},AxesOrigin->{-0.1,0},Frame->True,FrameLabel->{"\!\(\*SubsuperscriptBox[OverscriptBox[\(m\), \(^\)], \(PS\), \(2\)]\)",ylabel},LabelStyle->Directive[15,Black],PlotLegends->{Column[MaTeX[{"\\text{Extrapolation}", "\\text{(higher }\\beta)"},FontSize->20],Alignment->Left]}];
AppendTo[showplots,DashHighBeta];
LineHighBeta=Plot[{f[CombinedResultsHighBeta[[1,1]],CombinedResultsHighBeta[[2,1]],x,CombinedResultsHighBeta[[3,1]],0]+CombinedResultsHighBeta[[1,2]],f[CombinedResultsHighBeta[[1,1]],CombinedResultsHighBeta[[2,1]],x,CombinedResultsHighBeta[[3,1]],0]-CombinedResultsHighBeta[[1,2]]},{x,0,1.1Max[xvalues]},PlotStyle->{{Dashed,Gray,Opacity[0.0001]},{Dashed,Gray,Opacity[0.0001]}},Filling->True,AxesOrigin->{-0.1,0},Frame->True,FrameLabel->{"\!\(\*SubsuperscriptBox[OverscriptBox[\(m\), \(^\)], \(PS\), \(2\)]\)",ylabel},LabelStyle->Directive[22,Black]];
AppendTo[showplots,LineHighBeta];
Export[OutputFilenameBase<>"_highbeta.dat",CombinedResultsHighBeta,"Table","FieldSeparators"->" "],
Export[OutputFilenameBase<>"_highbeta.dat",Array["--"&,{4,2}],"Table","FieldSeparators"->" "]];


(* ::Input:: *)
FinalPlot=Show[showplots,PlotLabel->MaTeX["\\chi^2/\\text{d.o.f.}="<>ToString[NumberForm[ObservableFinal[[4]],{3,2}]],FontSize->24],Epilog->{Directive[Red,Opacity[0.1]],If[UseObservable,limit1,Nothing],If[UseObservable,limit2,Nothing]}];


(* ::Input:: *)
Export[OutputFilenameBase<>".pdf",FinalPlot,"PDF",AllowRasterization->False,ImageResolution->1200];


(* ::Input:: *)
If[UseObservable,{Export[OutputFilenameBase<>".dat",CombinedResultsFinal,"Table","FieldSeparators"->" "];},Export[OutputFilenameBase<>".dat",Array["--"&,{4,2}],"Table","FieldSeparators"->" "]];
