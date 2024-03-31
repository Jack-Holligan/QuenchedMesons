#!/usr/bin/env wolframscript
(* ::Package:: *)

SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1}];


(* ::Input:: *)
Needs["MaTeX`"]
Rep="_SED_REP_";
channel="_SED_CHANNEL_";
observable="_SED_OBSERVABLE_";
shortchannels=<| "pseudoscalar" -> "PS", "vector" -> "V", "axialvector" -> "AV", "scalar" -> "S", "tensor" -> "T", "axialtensor" -> "AT" |>;
shortchannel=shortchannels[channel];
Channels=<|"masses" -> <| "V" -> 1, "AV" -> 2, "S" -> 3, "T" -> 4, "AT" -> 5 |>, "decayconsts" -> <|"PS" -> 1, "V" -> 2, "AV" -> 3 |> |>;
Nc={"4","6","8"};
ChannelSymbol=If[Rep=="S","\\mathcal{"<>shortchannel<>"}","\\mathrm{"<>If[Rep=="AS",ToLowerCase[shortchannel],shortchannel]<>"}"];
Labels=If[observable=="decayconsts","\\frac{\\hat{f}^2_{\\rm "<>ChannelSymbol<>"}}{N_c^2}", "\\hat{m}^2_{"<>ChannelSymbol<>"}"];

fitforms=<|"F" -> a*x+b, "AS" -> a*x^2+b, "S" -> a*x^2+b |>;
fitform=fitforms[Rep];

(* ::Input:: *)
BaseDirectory=Directory[]<>"/processed_data";
OutputDirectory=BaseDirectory<>"/LargeN";


(* ::Section:: *)
(* Import data *)


(* ::Input:: *)
data=Table[(ArrayReshape[Import[BaseDirectory<>"/Sp"<>ToString[Nc[[n]]]<>"/Continuum/"<>Rep<>"/"<>Rep<>"_"<>observable<>".txt","Table"],{Length[Channels[observable]],2}]/.{x_,y_}->Around[x,y])[[Channels[observable][shortchannel]]],{n,1,Length[Nc]}];

(* ::Input:: *)
(*Import mass / decay data and square-root. Results are now Overscript[m, ^] / Overscript[f, ^].*)


(* ::Section:: *)
(* Computation *)


(* ::Input:: *)
Extrp=Select[Table[{1/ToExpression@Nc[[n]],data[[n]]},{n,1,Length[Nc]}],Length[#[[2]]]>0&];


(* ::Input:: *)
extrapdata=Map[{#[[1]],#[[2]]["Value"]}&,Extrp];


(* ::Input:: *)
extraperrors=Map[#[[2]]["Error"]&,Extrp];


(* ::Input:: *)
fit=NonlinearModelFit[extrapdata,fitform,{a,b},x,Weights->1/extraperrors^2,VarianceEstimatorFunction->(1&)];


(* ::Input:: *)
largeN=If[Length[Extrp]>2,Around[b/.fit["BestFitParameters"],fit["ParameterErrors"][[2]]],(Extrp[[1, 2]]/Extrp[[1, 1]]^2 - Extrp[[2, 2]]/Extrp[[2, 1]]^2)/(1/Extrp[[1, 1]]^2 - 1/Extrp[[2, 1]]^2)];


(* ::Input:: *)
ChiSq=If[Length[Extrp]>2,Sum[(fit[1/ToExpression[Nc[[i]]]]-Extrp[[i,2]]["Value"])^2/Extrp[[i,2]]["Error"]^2,{i,1,Length@Nc}],"--"];


(* ::Input:: *)
ExtrapolationPlot=Show[{ListPlot[Extrp,PlotStyle->Blue,PlotLegends->{MaTeX["\\text{Finite }N_c",FontSize->24]}],ListPlot[{{0, largeN}},PlotStyle->Black,PlotMarkers->"\[FilledSquare]",PlotLegends->MaTeX["\\text{Large-}N_c\\text{ limit}",FontSize->24]],Plot[fit[x],{x,0,0.3},PlotStyle->{Red,Dashed},PlotLegends->{MaTeX["\\text{Extrapolation}",FontSize->24]}]},Frame->True,LabelStyle->Directive[18,Black,FontFamily->Times],FrameLabel->{MaTeX["1/N_c",FontSize->24],MaTeX[Labels,FontSize->24]},RotateLabel->False,PlotLabel->MaTeX["\\chi^2/\\text{d.o.f.}="<>ToString[NumberForm[ChiSq,{3,2}]],FontSize->24],PlotRange->{{0,0.3},All},Axes->False];


(* ::Input:: *)
Export[OutputDirectory<>"/"<>Rep<>"_"<>channel<>"_"<>observable<>".pdf",ExtrapolationPlot,"PDF",AllowRasterization->False];


(* ::Section:: *)
(* Export *)


(* ::Input:: *)
Result={{shortchannel,largeN["Value"],largeN["Error"],ChiSq,If[Length[Extrp]>2,a/.fit["BestFitParameters"],"--"],If[Length[Extrp]>2,fit["ParameterErrors"][[1]],"--"]}};


(* ::Input:: *)
Export[OutputDirectory<>"/"<>Rep<>"_"<>channel<>"_"<>observable<>".txt",Result, "Table"];
