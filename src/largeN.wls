#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
Needs["MaTeX`"]
Rep="_SED_REP_";
channel="_SED_CHANNEL_";
observable="_SED_OBSERVABLE_";
shortchannels=<| "pseudoscalar" -> "PS", "vector" -> "V", "axialvector" -> "AV", "scalar" -> "S", "tensor" -> "T", "axialtensor" -> "AT" |>;
shortchannel=ChannelSymbols[channel];
Channels=<|"masses" -> <| "V" -> 1, "AV" -> 2, "S" -> 3, "T" -> 4, "AT" -> 5 |>, "decayconsts" -> <|"PS" -> 1, "V" -> 2, "AV" -> 3 |> |>;
Nc={"4","6","8"};
ChannelSymbol=If[Rep=="S","\\mathcal{"<>shortchannel<>"}","\\rm{"<>If[Rep=="AS",ToLowerCase[shortchannel],shortchannel]<>"}"];
Labels=If[observable=="decayconsts","\\frac{\\hat{f}^2_{\\rm "<>ChannelSymbol<>"}}{N_c^2}", "\\hat{m}^2_{"<>ChannelSymbol<>"}";


(* ::Input:: *)
BaseDirectory=Directory[]<>"/processed_data";
OutputDirectory=BaseDirectory<>"/LargeN";


(* ::Section:: *)
(* Import data *)


(* ::Input:: *)
data=(ArrayReshape[Import[BaseDirectory<>"/Sp"<>ToString[Nc[[n]]]<>"/Continuum/"<>Rep<>"/"<>Rep<>"_"<>observable<>".txt","Table"],{Length[Channels[observable]],2}]/.{x_,y_}->Around[x,y])/If[observable=="decayconsts",ToExpression[Nc[[n]]],1];


(* ::Input:: *)
(*Import mass / decay data and square-root. Results are now Overscript[m, ^] / Overscript[f, ^].*)


(* ::Section:: *)
(* Computation *)


(* ::Input:: *)
Extrp=Reap[Do[Sow[{1/ToExpression@Nc[[n]],data[[Channels[shortchannel]]]}],{n,1,Nc//Length}]][[2,1]];


(* ::Input:: *)
extrapdata=Extrp/.{x_,y_}->{x,y["Value"]};


(* ::Input:: *)
extraperrors=Extrp/.{x_,y_}->y["Error"];


(* ::Input:: *)
fit=NonlinearModelFit[extrapdata,a*x+b,{a,b},x,Weights->1/extraperrors^2,VarianceEstimatorFunction->(1&)];


(* ::Input:: *)
largeN=Around[b/.fit["BestFitParameters"],fitD["ParameterErrors"][[2]]];


(* ::Input:: *)
ChiSq=Sum[(fit[1/ToExpression[Nc[[i]]]]-Extrp[[i,2]]["Value"])^2/Extrp[[i,2]]["Error"]^2,{i,1,Length@Nc}];


(* ::Input:: *)
ExtrapolationPlot=Show[{ListPlot[Extrp,PlotStyle->Blue,PlotLegends->{MaTeX["\\text{Finite }N_c",FontSize->24]}],ListPlot[{{0, largeN}},PlotStyle->Black,PlotMarkers->"\[FilledSquare]",PlotLegends->MaTeX["\\text{Large-}N_c\\text{ limit}",FontSize->24]],Plot[fit[x],{x,0,0.3},PlotStyle->{Red,Dashed},PlotLegends->{MaTeX["\\text{Extrapolation}",FontSize->24]}]},Frame->True,LabelStyle->Directive[18,Black,FontFamily->Times],FrameLabel->{MaTeX["1/N_c",FontSize->24],MaTeX[Labels,FontSize->24]},RotateLabel->False,PlotLabel->MaTeX["\\chi^2/\\text{d.o.f.}="<>ToString[NumberForm[ChiSq,{3,2}]],FontSize->24],PlotRange->{{0,0.3},All},Axes->False];


(* ::Input:: *)
Export[OutputDirectory<>"/"<>Rep<>"_"<>channel<>"_"<>observable<>".pdf",ExtrapolationPlot,"PDF",AllowRasterization->False];


(* ::Section:: *)
(* Export *)


(* ::Input:: *)
Result={{shortchannel,largeN["Value"],largeN["Error"],ChiSq}};


(* ::Input:: *)
Export[OutputDirectory<>"/"<>Rep<>"_"<>channel"_"<>observable<>".txt",Result, "Table"];
