(* 
    -----------------------------------------------------------------
    -                         Dataset Import
    -----------------------------------------------------------------
*)

LaunchKernels[];
Print["The total number of cores is ", Length @ Kernels[]];

rawdata = Import["tracer_out.csv"];
initTime = rawdata[[1, 1]];
crtG=1.*^-15;
crtS=1.*^-8;
data = ReplacePart[rawdata[[#]], {
  3 -> Round[rawdata[[#, 3]]],
  1 -> (rawdata[[#, 1]] - initTime)
  (* 7 -> If[Abs[rawdata[[#, 7]]] < crtG, 0, rawdata[[#, 7]]],
  8 -> If[Abs[rawdata[[#, 8]]] < crtS, 0, rawdata[[#, 8]]] *)
  }]& /@ Range[Length @ rawdata];
time = DeleteDuplicates[data[[All, 1]]];
Print["The number of time points are ", Length @ time];
height = ReverseSort @ DeleteDuplicates[data[[All, 3]]];
Print["The heights are ", height];

(* 
    ------------------------------------------------------------------
    -                            Functions
    ------------------------------------------------------------------ 
*)

calTcG[time_, height_] :=
  Block[{dataCal, dataEs, solEs, REs, dataOz, solOz, ROz, posTime, posHeight,
     dataUg, solUg, RUg},
    posTime = Map[Boole[# == time]&, data[[All, 1]]];
    posHeight = Map[Boole[# == height]&, data[[All, 3]]];
    dataCal = Pick[data, posTime * posHeight, 1];
    (**)
    dataEs = {#2, #4}& @@@ dataCal;
    solEs = Interpolation[dataEs, InterpolationOrder -> 1] //
      SetPrecision[#, 16]&;
    (**)
    dataUg = {#2, #5}& @@@ dataCal;
    solUg = Interpolation[dataUg, InterpolationOrder -> 1] //
      SetPrecision[#, 16]&;
    (**)
    dataOz = {#2, #7}& @@@ dataCal;
    solOz = Interpolation[dataOz, InterpolationOrder -> 1] //
      SetPrecision[#, 16]&;
    NIntegrate[
      (1. - solEs[x]) * solUg[x] * solOz[x]
      , {x, 0., diameter}
      , WorkingPrecision -> 16
      , PrecisionGoal -> 10
      , AccuracyGoal -> 10
    ]
  ]

calTcS[time_, height_] :=
  Block[{dataCal, dataEs, solEs, REs, dataOz, solOz, ROz, posTime, posHeight,
     dataUs, solUs, RUg},
    posTime = Map[Boole[# == time]&, data[[All, 1]]];
    posHeight = Map[Boole[# == height]&, data[[All, 3]]];
    dataCal = Pick[data, posTime * posHeight, 1];
    (**)
    dataEs = {#2, #4}& @@@ dataCal;
    solEs = Interpolation[dataEs, InterpolationOrder -> 1] //
      SetPrecision[#, 16]&;
    (**)
    dataUs = {#2, #6}& @@@ dataCal;
    solUs = Interpolation[dataUs, InterpolationOrder -> 1] //
      SetPrecision[#, 16]&;
    (**)
    dataOz = {#2, #8}& @@@ dataCal;
    solOz = Interpolation[dataOz, InterpolationOrder -> 1] //
      SetPrecision[#, 16]&;
    NIntegrate[
      solEs[x] * solUs[x] * solOz[x]
      , {x, 0., diameter}
      , WorkingPrecision -> 16
      , PrecisionGoal -> 10
      , AccuracyGoal -> 10
    ]
  ]

(* 
    -----------------------------------------------------------------
    -                            Export
    -----------------------------------------------------------------
*)

diameter = 0.0762;

(* Oz0 = 4.0022943*^-8;
Tc0 = 6.6167027*^-5; *)

Print["Start calculating the data of gas tracer."];

Print[ Block[{}, tcG = ParallelTable[
  calTcG[time[[i]], height[[j]]]
  , {i, 1, Length @ time}, {j, 2, 3}
  , Method -> "CoarsestGrained"
   ]; ] // AbsoluteTiming ];

Print["Start calculating the data of solids tracer."];

Print[ Block[{}, tcS = ParallelTable[
  calTcS[time[[i]], height[[j]]]
  , {i, 1, Length @ time}, {j, 2, 3}
  , Method -> "CoarsestGrained"
   ]; ] // AbsoluteTiming ];

(* Export["data-time_series/tcG.csv", tcG, "CSV"];
Export["data-time_series/tcS.csv", tcS, "CSV"]; *)

dataExp = Transpose[{time, tcG[[All, 1]], tcG[[All, 2]], tcS[[All, 1]], tcS[[All, 2]]}];

Export["data-time_series/data-tc.csv", dataExp, "CSV"];

Print["The calculation is finished."];

(* 
  CMD line running code: 
    wolframscript -file '.\cal-tracer_time_series.wl' 
*)
