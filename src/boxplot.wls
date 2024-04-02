#!/usr/bin/env wolframscript
(* ::Package:: *)

SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1}];

Needs["MaTeX`"]
Basepath = Directory[] <> "/processed_data";
Outputpath = Basepath <> "/boxplots"
Nc = {4, 6, 8, "\\infty"};
Labels = {"Sp(4)", "Sp(6)", "Sp(8)", "Sp(\[Infinity])"};
Rep = {"F", "AS", "S"};
Representations = {"Fundamental", "Antisymmetric", "Symmetric"};
MassChannels = {"V", "AV", "S", "T", "AT"};
DecayChannels = {"PS", "V", "AV"};
Colours = {Purple, Darker[Green], Blue, Orange};
MassChannelsLaTeX = {MaTeX["\\rm V", FontSize -> 24],
   MaTeX[ "\\rm AV", FontSize -> 24], MaTeX[ "\\rm S", FontSize -> 24],
   MaTeX[ "\\rm T", FontSize -> 24], MaTeX["\\rm AT", FontSize -> 24]};
MassChannelsLaTeXAS = {MaTeX["\\rm v", FontSize -> 24],
   MaTeX[ "\\rm av", FontSize -> 24], MaTeX[ "\\rm s", FontSize -> 24],
   MaTeX[ "\\rm t", FontSize -> 24], MaTeX["\\rm at", FontSize -> 24]};
MassChannelsLaTeXS = {MaTeX["\\mathcal{V}", FontSize -> 24],
   MaTeX[ "\\mathcal{AV}", FontSize -> 24],
   MaTeX[ "\\mathcal{S}", FontSize -> 24],
   MaTeX[ "\\mathcal{T}", FontSize -> 24],
   MaTeX["\\mathcal{AT}", FontSize -> 24]};
DecayChannelsLaTeX = {MaTeX["\\rm PS", FontSize -> 24],
   MaTeX["\\rm V", FontSize -> 24], MaTeX[ "\\rm AV", FontSize -> 24]};
DecayChannelsLaTeXAS = {MaTeX["\\rm ps", FontSize -> 24],
   MaTeX["\\rm v", FontSize -> 24], MaTeX[ "\\rm av", FontSize -> 24]};
DecayChannelsLaTeXS = {MaTeX["\\mathcal{PS}", FontSize -> 24],
   MaTeX["\\mathcal{V}", FontSize -> 24],
   MaTeX[ "\\mathcal{AV}", FontSize -> 24]};

(* Decays *)

DecaysF =
  Reap[Do[Sow[
      ArrayReshape[
       Import[Basepath <> "/Sp" <> ToString[Nc[[i]]] <>
         "/continuum/F/F_decayconsts.txt",
        "Table"], {DecayChannels // Length, 2}]], {i, 1, 3}]][[2, 1]];
DecaysAS =
  Reap[Do[Sow[
      ArrayReshape[
       Import[Basepath <> "/Sp" <> ToString[Nc[[i]]] <>
         "/continuum/AS/AS_decayconsts.txt",
        "Table"], {DecayChannels // Length, 2}]], {i, 1, 3}]][[2, 1]];
DecaysS =
  Reap[Do[Sow[
      ArrayReshape[
       Import[Basepath <> "/Sp" <> ToString[Nc[[i]]] <>
         "/continuum/S/S_decayconsts.txt",
        "Table"], {DecayChannels // Length, 2}]], {i, 1, 3}]][[2, 1]];

Do[DecaysF[[n, ;;]]/Nc[[n]], {n, 1, 3}];
Do[DecaysAS[[n, ;;]]/Nc[[n]]^2, {n, 1, 3}]
Do[DecaysS[[n, ;;]]/Nc[[n]]^2, {n, 1, 3}]

DecaysNF =
  Import[Basepath <> "/largeN/F_decayconsts.txt",
    "Table"] /. {x_, y_, z_, w_, a_, da_} -> {Around[y, z]};
DecaysNAS =
  Import[Basepath <> "/largeN/AS_decayconsts.txt",
    "Table"] /. {x_, y_, z_, w_, a_, da_} -> {Around[y, z]};
DecaysNS =
  Import[Basepath <> "/largeN/S_decayconsts.txt",
    "Table"] /. {x_, y_, z_, w_, a_, da_} -> {Around[y, z]};

AppendTo[DecaysF, DecaysNF /. {x_} -> {x["Value"], x["Error"]}];
AppendTo[DecaysAS, DecaysNAS /. {x_} -> {x["Value"], x["Error"]}];
AppendTo[DecaysS, DecaysNS /. {x_} -> {x["Value"], x["Error"]}];

For[i = 1, i <= Dimensions[DecaysF][[1]], ++i,
 For[j = 1, j <= Dimensions[DecaysF][[2]], ++j,
  For[k = 1, k <= Dimensions[DecaysF][[3]], ++k,
   If[StringMatchQ[ToString[DecaysF[[i, j, k]]], "--"],
    DecaysF[[i, j]] = {100, 0}]]]]

DecaysFPlots =
  Reap[Do[Do[
      Sow[Plot[{DecaysF[[i, j, 1]] + DecaysF[[i, j, 2]],
         DecaysF[[i, j, 1]] -
          DecaysF[[i, j, 2]]}, {z, (Nc // Length)*(j - 1) + i -
          0.9, (Nc // Length)*(j - 1) + i - 0.1},
        Filling -> {1 -> {{2}, Colours[[i]]}},
        PlotStyle -> Colours[[i]]]], {j, 1, DecayChannels // Length}], {i, 1,
      Nc // Length}]][[2, 1]];

i = 1;
DecaysFPlots[[1]] =
  Plot[{DecaysF[[i, 1, 1]] + DecaysF[[i, 1, 2]],
    DecaysF[[i, 1, 1]] - DecaysF[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> {All, {-0.001, 0.1}}];
i = 2;
DecaysFPlots[[4]] =
  Plot[{DecaysF[[i, 1, 1]] + DecaysF[[i, 1, 2]],
    DecaysF[[i, 1, 1]] - DecaysF[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
i = 3;
DecaysFPlots[[7]] =
  Plot[{DecaysF[[i, 1, 1]] + DecaysF[[i, 1, 2]],
    DecaysF[[i, 1, 1]] - DecaysF[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
i = 4;
DecaysFPlots[[10]] =
  Plot[{DecaysF[[i, 1, 1]] + DecaysF[[i, 1, 2]],
    DecaysF[[i, 1, 1]] - DecaysF[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
AppendTo[DecaysFPlots,
  Table[Graphics[{Text[
      Style[DecayChannelsLaTeX[[i]], Black, FontFamily -> "Times",
       FontSlant -> Plain,
       FontSize -> 24], {(Nc // Length)*(i - 0.5), -0.005}]}], {i, 1,
    DecayChannelsLaTeX // Length}]];

FDecayPlot=Show[DecaysFPlots, PlotRange -> {All, {-0.01, 0.1}}, ImageSize -> 800,
 LabelStyle -> Directive[24, Black, FontFamily -> Times], Frame -> True,
 FrameLabel -> {"", MaTeX["\\frac{\\hat{f}^2_{\\chi}}{N_c}", FontSize -> 24]},
  RotateLabel -> False,
 GridLines -> {Table[
    i, {i, Nc // Length, (Nc // Length)*(DecayChannels // Length) - 1,
     Nc // Length}], Table[i, {i, 1/100, 0.5, 1/100}]},
 GridLinesStyle -> {Line, Dotted}, ImagePadding -> {{80, 10}, {3, None}},
 PlotLabel -> MaTeX["\\text{(f) fermions}", FontSize -> 24]]

Export[Outputpath <> "/FundamentalDecay.pdf", FDecayPlot, "PDF"];

For[i = 1, i <= Dimensions[DecaysAS][[1]], ++i,
 For[j = 1, j <= Dimensions[DecaysAS][[2]], ++j,
  For[k = 1, k <= Dimensions[DecaysAS][[3]], ++k,
   If[StringMatchQ[ToString[DecaysAS[[i, j, k]]], "--"],
    DecaysAS[[i, j]] = {100, 0}]]]]

DecaysASPlots =
  Reap[Do[Do[
      Sow[Plot[{DecaysAS[[i, j, 1]] + DecaysAS[[i, j, 2]],
         DecaysAS[[i, j, 1]] -
          DecaysAS[[i, j, 2]]}, {z, (Nc // Length)*(j - 1) + i -
          0.9, (Nc // Length)*(j - 1) + i - 0.1},
        Filling -> {1 -> {{2}, Colours[[i]]}},
        PlotStyle -> Colours[[i]]]], {j, 1, DecayChannels // Length}], {i, 1,
      Length[Nc]}]][[2, 1]];

i = 1;
DecaysASPlots[[1]] =
  Plot[{DecaysAS[[i, 1, 1]] + DecaysAS[[i, 1, 2]],
    DecaysAS[[i, 1, 1]] - DecaysAS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> {All, {-0.001, 0.1}}];
i = 2;
DecaysASPlots[[4]] =
  Plot[{DecaysAS[[i, 1, 1]] + DecaysAS[[i, 1, 2]],
    DecaysAS[[i, 1, 1]] - DecaysAS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> {All, {-0.001, 0.1}}];
i = 3;
DecaysASPlots[[7]] =
  Plot[{DecaysAS[[i, 1, 1]] + DecaysAS[[i, 1, 2]],
    DecaysAS[[i, 1, 1]] - DecaysAS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> {All, {-0.001, 0.1}}];
i = 4;
DecaysASPlots[[10]] =
  Plot[{DecaysAS[[i, 1, 1]] + DecaysAS[[i, 1, 2]],
    DecaysAS[[i, 1, 1]] - DecaysAS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
AppendTo[DecaysASPlots,
  Table[Graphics[{Text[
      Style[DecayChannelsLaTeXAS[[i]], Black, FontFamily -> "Times",
       FontSlant -> Plain,
       FontSize -> 24], {(Nc // Length)*(i - 0.5), -0.01}]}], {i, 1,
    DecayChannelsLaTeXAS // Length}]];

ASDecayPlot=Show[DecaysASPlots, PlotRange -> {All, {-0.02, 0.35}}, ImageSize -> 800,
 LabelStyle -> Directive[24, Black, FontFamily -> Times], Frame -> True,
 FrameLabel -> {"",
   MaTeX["\\frac{\\hat{f}^2_{\\chi}}{N_c^2}", FontSize -> 24]},
 RotateLabel -> False,
 GridLines -> {Table[
    i, {i, Nc // Length, (Nc // Length)*(DecayChannels // Length) - 1,
     Nc // Length}], Table[i, {i, 1/100, 0.5, 1/100}]},
 GridLinesStyle -> {Line, Dotted}, ImagePadding -> {{80, 10}, {3, None}},
 PlotLabel -> MaTeX["\\text{(as) fermions}", FontSize -> 24]]

Export[Outputpath <> "/AntisymmetricDecay.pdf", ASDecayPlot, "PDF"];

For[i = 1, i <= Dimensions[DecaysS][[1]], ++i,
 For[j = 1, j <= Dimensions[DecaysS][[2]], ++j,
  For[k = 1, k <= Dimensions[DecaysS][[3]], ++k,
   If[StringMatchQ[ToString[DecaysS[[i, j, k]]], "--"],
    DecaysS[[i, j]] = {100, 0}]]]]

DecaysSPlots =
  Reap[Do[Do[
      Sow[Plot[{DecaysS[[i, j, 1]] + DecaysS[[i, j, 2]],
         DecaysS[[i, j, 1]] -
          DecaysS[[i, j, 2]]}, {z, (Nc // Length)*(j - 1) + i -
          0.9, (Nc // Length)*(j - 1) + i - 0.1},
        Filling -> {1 -> {{2}, Colours[[i]]}},
        PlotStyle -> Colours[[i]]]], {j, 1, DecayChannels // Length}], {i, 1,
      Nc // Length}]][[2, 1]];

i = 1;
DecaysSPlots[[1]] =
  Plot[{DecaysS[[i, 1, 1]] + DecaysS[[i, 1, 2]],
    DecaysS[[i, 1, 1]] - DecaysS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> {All, {-0.01, 0.3}}];
i = 2;
DecaysSPlots[[4]] =
  Plot[{DecaysS[[i, 1, 1]] + DecaysS[[i, 1, 2]],
    DecaysS[[i, 1, 1]] - DecaysS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
i = 3;
DecaysSPlots[[7]] =
  Plot[{DecaysS[[i, 1, 1]] + DecaysS[[i, 1, 2]],
    DecaysS[[i, 1, 1]] - DecaysS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
i = 4;
DecaysSPlots[[10]] =
  Plot[{DecaysS[[i, 1, 1]] + DecaysS[[i, 1, 2]],
    DecaysS[[i, 1, 1]] - DecaysS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
AppendTo[DecaysSPlots,
  Table[Graphics[{Text[
      Style[DecayChannelsLaTeXS[[i]], Black, FontFamily -> "Times",
       FontSlant -> Plain,
       FontSize -> 24], {(Nc // Length)*(i - 0.5), -0.02}]}], {i, 1,
    DecayChannelsLaTeX // Length}]];

SDecayPlot=Show[DecaysSPlots, PlotRange -> {All, {-0.05, 0.8}}, ImageSize -> 800,
 LabelStyle -> Directive[24, Black, FontFamily -> Times], Frame -> True,
 FrameLabel -> {"",
   MaTeX["\\frac{\\hat{f}^2_{\\chi}}{N_c^2}", FontSize -> 24]},
 RotateLabel -> False,
 GridLines -> {Table[
    i, {i, Nc // Length, (Nc // Length)*(DecayChannels // Length),
     Nc // Length}], Table[i, {i, 0, 0.8, 1/20}]},
 GridLinesStyle -> {Line, Dotted}, ImagePadding -> {{80, 10}, {3, None}},
 PlotLabel -> MaTeX["\\text{(s) fermions}", FontSize -> 24]]

Export[Outputpath <> "/SymmetricDecay.pdf", SDecayPlot, "PDF"];

(* Masses *)

MassesF =
  Reap[Do[Sow[
      ArrayReshape[
       Import[Basepath <> "/Sp" <> ToString[Nc[[i]]] <>
         "/continuum/F/F_masses.txt",
        "Table"], {MassChannels // Length, 2}]], {i, 1, 3}]][[2, 1]];
MassesAS =
  Reap[Do[Sow[
      ArrayReshape[
       Import[Basepath <> "/Sp" <> ToString[Nc[[i]]] <>
         "/continuum/AS/AS_masses.txt",
        "Table"], {MassChannels // Length, 2}]], {i, 1, 3}]][[2, 1]];
MassesS =
  Reap[Do[Sow[
      ArrayReshape[
       Import[Basepath <> "/Sp" <> ToString[Nc[[i]]] <>
         "/continuum/S/S_masses.txt",
        "Table"], {MassChannels // Length, 2}]], {i, 1, 3}]][[2, 1]];

largeNF =
  Import[Basepath <> "/largeN/F_masses.txt",
    "Table"] /. {x_, y_, z_, w_, a_, da_} -> {Around[y, z]};
largeNAS =
  Import[Basepath <> "/largeN/AS_masses.txt",
    "Table"] /. {x_, y_, z_, w_, a_, da_} -> {Around[y, z]};
largeNS =
  Import[Basepath <> "/largeN/S_masses.txt",
    "Table"] /. {x_, y_, z_, w_, a_, da_} -> {Around[y, z]};

AppendTo[MassesF, largeNF /. {x_} -> {x["Value"], x["Error"]}];
AppendTo[MassesAS, largeNAS /. {x_} -> {x["Value"], x["Error"]}];
AppendTo[MassesS, largeNS /. {x_} -> {x["Value"], x["Error"]}];

For[i = 1, i <= Dimensions[MassesF][[1]], ++i,
 For[j = 1, j <= Dimensions[MassesF][[2]], ++j,
  For[k = 1, k <= Dimensions[MassesF][[3]], ++k,
   If[StringMatchQ[ToString[MassesF[[i, j, k]]], "--"],
    MassesF[[i, j]] = {100, 0}]]]]

MassesFPlots =
  Reap[Do[Do[
      Sow[Plot[{MassesF[[i, j, 1]] + MassesF[[i, j, 2]],
         MassesF[[i, j, 1]] -
          MassesF[[i, j, 2]]}, {z, (Nc // Length)*(j - 1) + i -
          0.9, (Nc // Length)*(j - 1) + i - 0.1},
        Filling -> {1 -> {{2}, Colours[[i]]}},
        PlotStyle -> Colours[[i]]]], {j, 1, MassChannels // Length}], {i, 1,
      Nc // Length}]][[2, 1]];

i = 1;
MassesFPlots[[1]] =
  Plot[{MassesF[[i, 1, 1]] + MassesF[[i, 1, 2]],
    MassesF[[i, 1, 1]] - MassesF[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> {All, {-0.001, 0.5}}];
i = 2;
MassesFPlots[[6]] =
  Plot[{MassesF[[i, 1, 1]] + MassesF[[i, 1, 2]],
    MassesF[[i, 1, 1]] - MassesF[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
i = 3;
MassesFPlots[[11]] =
  Plot[{MassesF[[i, 1, 1]] + MassesF[[i, 1, 2]],
    MassesF[[i, 1, 1]] - MassesF[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
i = 4;
MassesFPlots[[16]] =
  Plot[{MassesF[[i, 1, 1]] + MassesF[[i, 1, 2]],
    MassesF[[i, 1, 1]] - MassesF[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
AppendTo[MassesFPlots,
  Table[Graphics[{Text[
      Style[MassChannelsLaTeX[[i]], Black, FontFamily -> "Times",
       FontSlant -> Plain,
       FontSize -> 24], {(Nc // Length)*(i - 0.5), -0.1}]}], {i, 1,
    MassChannelsLaTeX // Length}]];

FMassPlot=Show[MassesFPlots, PlotRange -> {All, {-0.2, 1.5}}, ImageSize -> 800,
 LabelStyle -> Directive[24, Black, FontFamily -> Times], Frame -> True,
 FrameLabel -> {"", MaTeX["\\hat{m}^2_{\\chi}", FontSize -> 24]},
 RotateLabel -> False,
 GridLines -> {Table[
    i, {i, Nc // Length, (Nc // Length)*(MassChannels // Length) - 1,
     Nc // Length}], Table[i, {i, 1/10, 1.5, 1/10}]},
 GridLinesStyle -> {Line, Dotted}, ImagePadding -> {{80, 10}, {3, None}},
 PlotLabel -> MaTeX["\\text{(f) fermions}", FontSize -> 24]]

Export[Outputpath <> "/FundamentalMass.pdf", FMassPlot, "PDF"];

For[i = 1, i <= Dimensions[MassesAS][[1]], ++i,
 For[j = 1, j <= Dimensions[MassesAS][[2]], ++j,
  For[k = 1, k <= Dimensions[MassesAS][[3]], ++k,
   If[StringMatchQ[ToString[MassesAS[[i, j, k]]], "--"],
    MassesAS[[i, j]] = {100, 0}]]]]

MassesASPlots =
  Reap[Do[Do[
      Sow[Plot[{MassesAS[[i, j, 1]] + MassesAS[[i, j, 2]],
         MassesAS[[i, j, 1]] -
          MassesAS[[i, j, 2]]}, {z, (Nc // Length)*(j - 1) + i -
          0.9, (Nc // Length)*(j - 1) + i - 0.1},
        Filling -> {1 -> {{2}, Colours[[i]]}},
        PlotStyle -> Colours[[i]]]], {j, 1, MassChannels // Length}], {i, 1,
      Nc // Length}]][[2, 1]];

i = 1;
MassesASPlots[[1]] =
  Plot[{MassesAS[[i, 1, 1]] + MassesAS[[i, 1, 2]],
    MassesAS[[i, 1, 1]] - MassesAS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> {All, {-0.001, 2.5}}];
i = 2;
MassesASPlots[[6]] =
  Plot[{MassesAS[[i, 1, 1]] + MassesAS[[i, 1, 2]],
    MassesAS[[i, 1, 1]] - MassesAS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
i = 3;
MassesASPlots[[11]] =
  Plot[{MassesAS[[i, 1, 1]] + MassesAS[[i, 1, 2]],
    MassesAS[[i, 1, 1]] - MassesAS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
i = 4;
MassesASPlots[[16]] =
  Plot[{MassesAS[[i, 1, 1]] + MassesAS[[i, 1, 2]],
    MassesAS[[i, 1, 1]] - MassesAS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
AppendTo[MassesASPlots,
  Table[Graphics[{Text[
      Style[MassChannelsLaTeXAS[[i]], Black, FontFamily -> "Times",
       FontSlant -> Plain,
       FontSize -> 24], {(Nc // Length)*(i - 0.5), -0.15}]}], {i, 1,
    MassChannelsLaTeX // Length}]];

ASMassPlot=Show[MassesASPlots, PlotRange -> {All, {-0.3, 2.5}}, ImageSize -> 800,
 LabelStyle -> Directive[24, Black, FontFamily -> Times], Frame -> True,
 FrameLabel -> {"", MaTeX["\\hat{m}^2_{\\chi}", FontSize -> 24]},
 RotateLabel -> False,
 GridLines -> {Table[
    i, {i, Nc // Length, (Nc // Length)*(MassChannels // Length),
     Nc // Length}], Table[i, {i, 1/10, 2.5, 1/10}]},
 GridLinesStyle -> {Line, Dotted}, ImagePadding -> {{80, 10}, {3, None}},
 PlotLabel -> MaTeX["\\text{(as) fermions}", FontSize -> 24]]

Export[Outputpath <> "/AntisymmetricMass.pdf", ASMassPlot, "PDF"];

For[i = 1, i <= Dimensions[MassesS][[1]], ++i,
 For[j = 1, j <= Dimensions[MassesS][[2]], ++j,
  For[k = 1, k <= Dimensions[MassesS][[3]], ++k,
   If[StringMatchQ[ToString[MassesS[[i, j, k]]], "--"],
    MassesS[[i, j]] = {100, 0}]]]]

MassesSPlots =
  Reap[Do[Do[
      Sow[Plot[{MassesS[[i, j, 1]] + MassesS[[i, j, 2]],
         MassesS[[i, j, 1]] -
          MassesS[[i, j, 2]]}, {z, (Nc // Length)*(j - 1) + i -
          0.9, (Nc // Length)*(j - 1) + i - 0.1},
        Filling -> {1 -> {{2}, Colours[[i]]}},
        PlotStyle -> Colours[[i]]]], {j, 1, MassChannels // Length}], {i, 1,
      Nc // Length}]][[2, 1]];

i = 1;
MassesSPlots[[1]] =
  Plot[{MassesS[[i, 1, 1]] + MassesS[[i, 1, 2]],
    MassesS[[i, 1, 1]] - MassesS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> {All, {-0.001, 2.5}}];
i = 2;
MassesSPlots[[6]] =
  Plot[{MassesS[[i, 1, 1]] + MassesS[[i, 1, 2]],
    MassesS[[i, 1, 1]] - MassesS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
i = 3;
MassesSPlots[[11]] =
  Plot[{MassesS[[i, 1, 1]] + MassesS[[i, 1, 2]],
    MassesS[[i, 1, 1]] - MassesS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
AppendTo[MassesSPlots,
  Table[Graphics[{Text[
      Style[MassChannelsLaTeXS[[i]], Black, FontFamily -> "Times",
       FontSlant -> Plain,
       FontSize -> 24], {(Nc // Length)*(i - 0.5), -0.3}]}], {i, 1,
    MassChannelsLaTeX // Length}]];
i = 4;
MassesSPlots[[16]] =
  Plot[{MassesS[[i, 1, 1]] + MassesS[[i, 1, 2]],
    MassesS[[i, 1, 1]] - MassesS[[i, 1, 2]]}, {z, i - 0.9, i - 0.1},
   Filling -> {1 -> {{2}, Colours[[i]]}}, PlotStyle -> Colours[[i]],
   PlotLegends -> {MaTeX["\\rm Sp(" <> ToString[Nc[[i]]] <> ")",
      FontSize -> 24]}, PlotRange -> Automatic];
AppendTo[MassesSPlots,
  Table[Graphics[{Text[
      Style[MassChannelsLaTeXS[[i]], Black, FontFamily -> "Times",
       FontSlant -> Plain,
       FontSize -> 24], {(Nc // Length)*(i - 0.5), -0.3}]}], {i, 1,
    MassChannelsLaTeX // Length}]];

SMassPlot=Show[MassesSPlots, PlotRange -> {{0, 19.5}, {-0.6, 5}}, ImageSize -> 800,
 LabelStyle -> Directive[24, Black, FontFamily -> Times], Frame -> True,
 FrameLabel -> {"", MaTeX["\\hat{m}^2_{\\chi}", FontSize -> 24]},
 RotateLabel -> False,
 GridLines -> {Table[
    i, {i, Nc // Length, (Nc // Length)*(MassChannels // Length),
     Nc // Length}], Table[i, {i, 1/5, 5, 1/5}]},
 GridLinesStyle -> {Line, Dotted}, ImagePadding -> {{80, 10}, {3, None}},
 PlotLabel -> MaTeX["\\text{(s) fermions}", FontSize -> 24]]

Export[Outputpath <> "/SymmetricMass.pdf", SMassPlot, "PDF"];

data1PS = DecaysF[[;; , 1, ;;]] /. {x_, y_} -> Around[x, y];
data1V = DecaysF[[;; , 2, ;;]] /. {x_, y_} -> Around[x, y];
data1AV = DecaysF[[;; , 3, ;;]] /. {x_, y_} -> Around[x, y];
data2 = MassesF[[;; , 1, ;;]] /. {x_, y_} -> Around[x, y];
RatioPS = Transpose@{1/Nc, (data1PS/data2)^-2};
RatioV = Transpose@{1/Nc, (data1V/data2)^-2};
RatioAV = Transpose@{1/Nc, (data1AV/data2)^-2};

fromFplot=ListPlot[{RatioPS, RatioV, RatioAV}, ImageSize -> 800, PlotRange -> All,
 PlotStyle -> Colours, Frame -> True,
 LabelStyle -> Directive[24, Black, FontFamily -> Times],
 FrameLabel -> {MaTeX["1/N_c", FontSize -> 24],
   MaTeX["\\frac{m^2_{\\rm V}}{f^2}", FontSize -> 24]}, RotateLabel -> False,
 PlotLegends -> {MaTeX["f_{\\text{PS}}", FontSize -> 24],
   MaTeX["f_{\\text{V}}", FontSize -> 24],
   MaTeX["f_{\\text{AV}}", FontSize -> 24]},
 PlotLabel -> MaTeX["\\text{Fundamental Representation}", FontSize -> 24]]

Export[Outputpath <> "/ftomRatioF.pdf", ftomFplot, "PDF"];

data1PS = DecaysAS[[;; , 1, ;;]] /. {x_, y_} -> Around[x, y];
data1V = DecaysAS[[;; , 2, ;;]] /. {x_, y_} -> Around[x, y];
data1AV = DecaysAS[[;; , 3, ;;]] /. {x_, y_} -> Around[x, y];
data2 = MassesAS[[;; , 1, ;;]] /. {x_, y_} -> Around[x, y];
RatioPS = Transpose@{1/Nc, (data1PS/data2)^-2};
RatioV = Transpose@{1/Nc, (data1V/data2)^-2};
RatioAV = Transpose@{1/Nc, (data1AV/data2)^-2};

fromASplot=ListPlot[{RatioPS, RatioV, RatioAV}, ImageSize -> 800, PlotRange -> All,
 PlotStyle -> Colours, Frame -> True,
 LabelStyle -> Directive[24, Black, FontFamily -> Times],
 FrameLabel -> {MaTeX["1/N_c", FontSize -> 24],
   MaTeX["\\frac{m^2_{\\rm v}}{f^2}", FontSize -> 24]}, RotateLabel -> False,
 PlotLegends -> {MaTeX["f_{\\text{ps}}", FontSize -> 24],
   MaTeX["f_{\\text{v}}", FontSize -> 24],
   MaTeX["f_{\\text{av}}", FontSize -> 24]},
 PlotLabel -> MaTeX["\\text{Antisymmetric Representation}", FontSize -> 24]]

Export[Outputpath <> "/ftomRatioAS.pdf", ftomASplot, "PDF"];

data1PS = DecaysS[[;; , 1, ;;]] /. {x_, y_} -> Around[x, y];
data1V = DecaysS[[;; , 2, ;;]] /. {x_, y_} -> Around[x, y];
data1AV = DecaysS[[;; , 3, ;;]] /. {x_, y_} -> Around[x, y];
data2 = MassesS[[;; , 1, ;;]] /. {x_, y_} -> Around[x, y];
RatioPS = Transpose@{1/Nc, (data1PS/data2)^-2};
RatioV = Transpose@{1/Nc, (data1V/data2)^-2};
RatioAV = Transpose@{1/Nc, (data1AV/data2)^-2};

ftomSplot=ListPlot[{RatioPS, RatioV, RatioAV}, ImageSize -> 800, PlotRange -> All,
 PlotStyle -> Colours, Frame -> True,
 LabelStyle -> Directive[24, Black, FontFamily -> Times],
 FrameLabel -> {MaTeX["1/N_c", FontSize -> 24],
   MaTeX["\\frac{m^2_{\\mathcal{V}}}{f^2}", FontSize -> 24]},
 RotateLabel -> False,
 PlotLegends -> {MaTeX["f_{\\mathcal{PS}}", FontSize -> 24],
   MaTeX["f_{\\mathcal{V}}", FontSize -> 24],
   MaTeX["f_{\\mathcal{AV}}", FontSize -> 24]},
 PlotLabel -> MaTeX["\\text{Symmetric Representation}", FontSize -> 24]]

Export[Outputpath <> "/ftomRatioS.pdf", ftomSplot, "PDF"];
