#!/usr/bin/env wolframscript
(* ::Package:: *)

Nc := _SED_NC_;
BaseDir := Directory[] <> "/processed_data/Sp" <> ToString[Nc] <> "/continuum";
OutputDir := BaseDir;

(* ::Input:: *)
Needs["MaTeX`"]
Reps={"F", "AS", "S"};
Representations={"Fundamental", "Antisymmetric", "Symmetric"};
MassChannels={"V", "AV", "S", "T", "AT"};
MassChannelsLaTeX={MaTeX["\\rm V",FontSize->20],MaTeX[ "\\rm AV",FontSize->20],MaTeX[ "\\rm S",FontSize->20],MaTeX[ "\\rm T",FontSize->20], MaTeX["\\rm AT",FontSize->20]};
DecayChannels={"PS","V", "AV"};
DecayChannelsLaTeX={MaTeX["\\rm PS",FontSize->20],MaTeX["\\rm V",FontSize->20],MaTeX[ "\\rm AV",FontSize->20]};
Colours={Black,Red,Darker[Green]};
MassPlots={};Masses={};Decays={};DecayPlots={};
For[i=1,i<=Length[Representations],++i,
AppendTo[Masses,ArrayReshape[Import[StringJoin[BaseDir,"/",Reps[[i]],"/",Reps[[i]],"_masses.txt"],"Table"],{Length[MassChannels],2}]]]
For[i=1,i<=Length[Representations],++i,
AppendTo[Decays,ArrayReshape[Import[StringJoin[BaseDir,"/",Reps[[i]],"/",Reps[[i]],"_decayconsts.txt"],"Table"],{Length[DecayChannels],2}]]]
(*Element={Representation, Channel, Value/Error}*)
For[i=1,i<=Dimensions[Masses][[1]],++i,For[j=1,j<=Dimensions[Masses][[2]],++j,For[k=1,k<=Dimensions[Masses][[3]],++k,If[StringMatchQ[ToString[Masses[[i,j,k]]],"--"],Masses[[i,j,k]]=0.]]]]
For[i=1,i<=Dimensions[Decays][[1]],++i,For[j=1,j<=Dimensions[Decays][[2]],++j,For[k=1,k<=Dimensions[Decays][[3]],++k,If[StringMatchQ[ToString[Decays[[i,j,k]]],"--"],Decays[[i,j,k]]=0.]]]]
MaxMass=Max[Masses[[;;,;;,1]]+Masses[[;;,;;,2]]];
MaxDecay=Max[Decays[[;;,;;,1]]+Decays[[;;,;;,2]]];


(* ::Input:: *)
Massplots=Reap[Do[Do[Sow[Plot[{Masses[[i,j,1]]+Masses[[i,j,2]],Masses[[i,j,1]]-Masses[[i,j,2]]},{z,3*(j-1)+i-0.9,3*(j-1)+i-0.1},Filling->{1->{{2},Colours[[i]]}},PlotStyle->Colours[[i]]]],{i,1,3}],{j,1,5}]][[2,1]];
i=1;Massplots[[i]]=Plot[{Masses[[i,1,1]]+Masses[[i,1,2]],Masses[[i,1,1]]-Masses[[i,1,2]]},{z,i-0.9,i-0.1},Filling->{1->{{2},Colours[[i]]}},PlotStyle->Colours[[i]],PlotLegends->{MaTeX["\\rm Fundamental",FontSize->18]},PlotRange->{All,{-0.2,5}}];
i=2;Massplots[[i]]=Plot[{Masses[[i,1,1]]+Masses[[i,1,2]],Masses[[i,1,1]]-Masses[[i,1,2]]},{z,i-0.9,i-0.1},Filling->{1->{{2},Colours[[i]]}},PlotStyle->Colours[[i]],PlotLegends->{MaTeX["\\rm Antiymmetric",FontSize->18]}];
i=3;Massplots[[i]]=Plot[{Masses[[i,1,1]]+Masses[[i,1,2]],Masses[[i,1,1]]-Masses[[i,1,2]]},{z,i-0.9,i-0.1},Filling->{1->{{2},Colours[[i]]}},PlotStyle->Colours[[i]],PlotLegends->{MaTeX["\\rm Symmetric",FontSize->18]}];
AppendTo[Massplots,Table[ Graphics[{Text[Style[MassChannelsLaTeX[[i]],Black,FontFamily->"Times",FontSlant->Plain,FontSize->20],{3*(i-1)+1.5,0}]}],{i,1,Dimensions[Masses][[2]]}]];
FinalMassplots=Show[Massplots,PlotRange->{All,{-0.2,5}},LabelStyle->Directive[18,Black,FontFamily->Times],Frame->True,FrameLabel->{"",MaTeX["\\hat{m}^2_{\\chi}",FontSize->18]},RotateLabel->False,GridLines->{Table[i,{i,3,14,3}],Table[i,{i,1/3,5,1/3}]},GridLinesStyle->{Line,Dotted},ImagePadding->{{80,10},{3,None}},PlotLabel->MaTeX["\\text{Sp("<>ToString[Nc]<>") chiral meson masses}",FontSize->20]]


(* ::Input:: *)
Export[StringJoin[OutputDir,"/chiral_mass_Sp",ToString[Nc],".pdf"],FinalMassplots,"PDF",AllowRasterization->False];


(* ::Input:: *)
Decayplots=Reap[Do[Do[Sow[Plot[{Decays[[i,j,1]]+Decays[[i,j,2]],Decays[[i,j,1]]-Decays[[i,j,2]]},{z,3*(j-1)+i-0.9,3*(j-1)+i-0.1},Filling->{1->{{2},Colours[[i]]}},PlotStyle->Colours[[i]]]],{i,1,3}],{j,1,3}]][[2,1]];
i=1;Decayplots[[i]]=Plot[{Decays[[i,1,1]]+Decays[[i,1,2]],Decays[[i,1,1]]-Decays[[i,1,2]]},{z,i-0.9,i-0.1},Filling->{1->{{2},Colours[[i]]}},PlotStyle->Colours[[i]],PlotLegends->{MaTeX["\\rm Fundamental",FontSize->18]},PlotRange->{All,{-0.01,0.2}}];
i=2;Decayplots[[i]]=Plot[{Decays[[i,1,1]]+Decays[[i,1,2]],Decays[[i,1,1]]-Decays[[i,1,2]]},{z,i-0.9,i-0.1},Filling->{1->{{2},Colours[[i]]}},PlotStyle->Colours[[i]],PlotLegends->{MaTeX["\\rm Antisymmetric",FontSize->18]}];
i=3;Decayplots[[i]]=Plot[{Decays[[i,1,1]]+Decays[[i,1,2]],Decays[[i,1,1]]-Decays[[i,1,2]]},{z,i-0.9,i-0.1},Filling->{1->{{2},Colours[[i]]}},PlotStyle->Colours[[i]],PlotLegends->{MaTeX["\\rm Symmetric",FontSize->18]}];
AppendTo[Decayplots,Table[Graphics[{Text[Style[DecayChannelsLaTeX[[i]],Black,FontFamily->"Times",FontSlant->Plain,FontSize->20],{3*(i-1)+1.5,-0.005}]}],{i,1,Dimensions[Decays][[2]]}]];
FinalDecayplots=Show[Decayplots,PlotRange->{All,{-0.01,0.2}},LabelStyle->Directive[18,Black,FontFamily->Times],Frame->True,FrameLabel->{"",MaTeX["\\hat{f}^2_{\\chi}",FontSize->18]},RotateLabel->False,GridLines->{Table[i,{i,3,8,3}],Table[i,{i,1/60,0.2,1/60}]},GridLinesStyle->{Line,Dotted},ImagePadding->{{80,10},{3,None}},PlotLabel->MaTeX["\\text{Sp("<>ToString[Nc]<>") chiral meson decay constants}",FontSize->20]]


(* ::Input:: *)
Export[StringJoin[NotebookDirectory[],"chiral_decayconst_Sp"<>ToString[Nc]<>".pdf"],FinalDecayplots,"PDF",AllowRasterization->False];
