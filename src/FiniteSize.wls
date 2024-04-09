#!/usr/bin/env wolframscript
(* ::Package:: *)

SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1}];

(* ::Input:: *)
Needs["MaTeX`"]
Basepath=Directory[]<>"/_SED_BASEPATH_";
Outputpath=Directory[]<>"/_SED_OUTPUTPATH_";
Directories={"S12T24B15.6","S16T32B15.6", "S20T40B15.6","S24T48B15.6"};
Volumes={12,16,20,24};
baremasses={"-0.8","-0.81","-0.82"};
Colours={Red,Blue,Darker[Green],Orange,Purple};
Rep="F";
Representation="Fundamental";
Allmasses={};
For[mass=1,mass<=Length[baremasses],++mass,{onemass={};
For[dir=1,dir<=Length[Directories],++dir,data=Import[Basepath<>"/"<>Directories[[dir]]<>"_m"<>Rep<>baremasses[[mass]]<>"/"<>"output_pseudoscalar.txt","Table"];
AppendTo[onemass,{NumberForm[baremasses[[mass]]],Volumes[[dir]],data[[3,1]],data[[4,1]]}/.{a_,b_,c_,d_}->{a,b,Sqrt[Around[c,d]]}]]};
AppendTo[Allmasses,onemass];]
FiniteSize[minf_,A_,L_]:=minf(1+A Exp[-minf*L]/(minf*L)^(3/2))


(* ::Input:: *)
Plotter[mass_]:=Module[{a1=mass,tmpdata,lplot,pplot,mplot,minf,nlm,scaledData},tmpdata=Allmasses[[a1,;;,2;;]];
nlm=NonlinearModelFit[tmpdata/.{x_,y_}->{x,y["Value"]},FiniteSize[m,A,L],{m,A},L,Weights->tmpdata/.{x_,y_}->y["Error"]^(-2),VarianceEstimatorFunction->(1&)];
minf=Around[m/.nlm["BestFitParameters"],nlm["ParameterErrors"][[1]]];
scaledData=tmpdata/.{x_,y_}->{x*minf["Value"],y/minf["Value"]};
lplot=ListPlot[scaledData,PlotStyle->Colours[[a1]],LabelStyle->Directive[18,Black,FontFamily->Times]];
pplot=Plot[nlm[x/minf["Value"]]/minf["Value"],{x,0,Max[Volumes]},PlotStyle->{Colours[[a1]],Dashed}];
Show[lplot,pplot,ImageSize->800,Frame->True,FrameLabel->{MaTeX["m_{\\infty}L",FontSize->30],MaTeX["\\frac{m_{\\text{PS}}}{m_{\\infty}}",FontSize->30]},PlotLabel->MaTeX["am_0 = "<>baremasses[[a1]]<>".~\\beta = 15.6",FontSize->30]]]


(* ::Input:: *)
Do[Export[Outputpath<>"/Sp6_beta15.6_m"<>Rep<>ToString[baremasses[[i]]]<>".pdf",Plotter[i],"PDF",AllowRasterization->False],{i,3}];
