#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*Needs["MaTeX`"]*)
(*Rep="_SED_REP_";*)
(*channel=_SED _CHANNEL _;*)
(*observable="_SED_OBSERVABLE_";*)
(**)
(*MassChannels={"V", "AV", "S", "T", "AT"};*)
(*massesLabelsF={"\\hat{m}^2_{\\rm V}","\\hat{m}^2_{\\rm AV}","\\hat{m}^2_{\\rm S}","\\hat{m}^2_{\\rm T}","\\hat{m}^2_{\\rm AT}"};*)
(*massesLabelsAS={"\\hat{m}^2_{\\rm v}","\\hat{m}^2_{\\rm av}","\\hat{m}^2_{\\rm s}","\\hat{m}^2_{\\rm t}","\\hat{m}^2_{\\rm at}"};*)
(*massesLabelsS={"\\hat{m}^2_{\\mathcal{V}}","\\hat{m}^2_{\\mathcal{AV}}","\\hat{m}^2_{\\mathcal{S}}","\\hat{m}^2_{\\mathcal{T}}","\\hat{m}^2_{\\mathcal{AT}}"};*)
(*DecayChannels={"PS","V", "AV"};*)
(*decayconstsLabelsF={"\\frac{\\hat{f}^2_{\\rm PS}}{N_c}","\\frac{\\hat{f}^2_{\\rm V}}{N_c}","\\frac{\\hat{f}^2_{\\rm AV}}{N_c}"};*)
(*decayconstsLabelsAS={"\\frac{\\hat{f}^2_{\\rm ps}}{N_c^2}","\\frac{\\hat{f}^2_{\\rm v}}{N_c^2}","\\frac{\\hat{f}^2_{\\rm av}}{N_c^2}"};*)
(*decayconstsLabelsS={"\\frac{\\hat{f}^2_{\\mathcal{PS}}}{N_c^2}","\\frac{\\hat{f}^2_{\\mathcal{V}}}{N_c^2}","\\frac{\\hat{f}^2_{\\mathcal{AV}}}{N_c^2}"};*)
(*Nc={"4","6","8"};*)


(* ::Input:: *)
(*BaseDirectory=Directory[]<>"/processed_data";*)
(*OutputDirectory=BaseDirectory<>"/LargeN";*)


(* ::Input:: *)
(*Channels=If[observable=="decayconsts",DecayChannels,MassChannels];*)


(* ::Section:: *)
(*Import data*)


(* ::Input:: *)
(*data=(ArrayReshape[Import[BaseDirectory<>"/sp"<>ToString[Nc[[n]]]<>"/Continuum/"<>Rep<>"/"<>Rep<>"_decayconsts.txt","Table"],{Channels//Length,2}]/.{x_,y_}->Around[x,y])/If[observable=="decayconsts",ToExpression[Nc[[n]]],1];*)


(* ::Input:: *)
(*(*Import mass / decay data and square-root. Results are now Overscript[m, ^] / Overscript[f, ^].*)*)


(* ::Section:: *)
(*Computation*)


(* ::Input:: *)
(*Extrp=Reap[Do[Sow[{1/ToExpression@Nc[[n]],data[[channel]]}],{n,1,Nc//Length}]][[2,1]];*)


(* ::Input:: *)
(*extrapdata=Extrp/.{x_,y_}->{x,y["Value"]};*)


(* ::Input:: *)
(*extraperrors=Extrp/.{x_,y_}->y["Error"];*)


(* ::Input:: *)
(*fit=NonlinearModelFit[extrapdata,a*x+b,{a,b},x,Weights->1/extraperrors^2,VarianceEstimatorFunction->(1&)];*)


(* ::Input:: *)
(*largeN=Around[b/.fit["BestFitParameters"],fitD["ParameterErrors"][[2]]];*)


(* ::Input:: *)
(*ChiSq=Sum[(fit[1/ToExpression[Nc[[i]]]]-Extrp[[i,2]]["Value"])^2/Extrp[[i,2]]["Error"]^2,{i,1,Length@Nc}];*)


(* ::Input:: *)
(*Labels=Evaluate[Symbol[observable<>"Labels"<>Rep]];*)


(* ::Input:: *)
(*ExtrapolationPlot=Show[*)
(*{*)
(*ListPlot[Extrp,PlotStyle->Blue,PlotLegends->{MaTeX["\\text{Finite }N_c",FontSize->24]}],*)
(*ListPlot[{{0, largeN}},PlotStyle->Black,PlotMarkers->"\[FilledSquare]",PlotLegends->MaTeX["\\text{Large-}N_c\\text{ limit}",FontSize->24]],*)
(*Plot[fit[x],{x,0,0.3},PlotStyle->{Red,Dashed},PlotLegends->{MaTeX["\\text{Extrapolation}",FontSize->24]}]*)
(*},*)
(*Frame->True,*)
(*LabelStyle->Directive[18,Black,FontFamily->Times],*)
(*FrameLabel->{MaTeX["1/N_c",FontSize->24],MaTeX[Labels,FontSize->24]},*)
(*RotateLabel->False,*)
(*PlotLabel->MaTeX["\\chi^2/\\text{d.o.f.}="<>ToString[NumberForm[ChiSq,{3,2}]],FontSize->24],*)
(*PlotRange->{{0,0.3},All},*)
(*Axes->False*)
(*];*)


(* ::Input:: *)
(*Export[OutputDirectory<>"/"<>Rep<>"_"<>ToString[DecayChannels[[channel]]]<>"_decay.PDF",ExtrapolationPlot,"PDF",AllowRasterization->False];*)


(* ::Section:: *)
(*Export*)


(* ::Input:: *)
(*Result={{Channels[[channel]],largeN["Value"],largeN["Error"],ChiSq}};*)


(* ::Input:: *)
(*Export[OutputDirectory<>"/"<>Rep<>"_"<>Channels[[channel]]"_"<>observable<>".txt",Result, "Table"];*)
