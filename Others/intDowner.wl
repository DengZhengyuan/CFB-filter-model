(* ::Package:: *)

BeginPackage["intDowner`"];


smooth::usage = "Function to smooth the profiles.";
findPosition::usage = "Find data positions that share the same height.";
calREs::usage = "Calculate the cross-sectional averaged solids holdup at one height and one time point.";
calRVs::usage = "Calculate the cross-sectional averaged solids velocity at one height and one time point";
calROz::usage = "Calculate the cross-sectional averaged ozone concentration at one height and one time point";


Begin["`Private`"];


smooth[data_, meanNum_] :=
    Module[{prcData, fltData},
        prcData = data[[All, 1]];
        fltData = MeanFilter[prcData, meanNum];
        Transpose[{fltData, #2& @@@ data}]
    ]


findPosition[data_, yC_] :=
    Block[{roundY},
        roundY = Round[data[[All, 2]], .01];
        Boole[# == yC]& /@ roundY
    ]


calREs[onedata_, yCoord_, diameter_] :=
    Module[{sol, calData, dataGrouped, res, pos},
        calData = onedata[[All, {1, 8}]];
        pos = findPosition[onedata,#]& /@ yCoord;
        dataGrouped = Pick[calData, #, 1]& /@ pos;
        dataGrouped = DeleteDuplicatesBy[#, First]& /@ dataGrouped;
        sol = Map[Interpolation[#, InterpolationOrder -> 1]&, dataGrouped];
        res = 1. / diameter * Quiet @ NIntegrate[#[x], {x, 0., diameter}]& /@ sol;
        res
    ]


calRVs[onedata_, yCoord_, diameter_] :=
    Module[{sol, calData, dataGrouped, pos, res},
        calData = {#1, #7 * #8}& @@@ onedata;
        pos = findPosition[onedata,#]& /@ yCoord;
        dataGrouped = Pick[calData, #, 1]& /@ pos;
        dataGrouped = DeleteDuplicatesBy[#, First]& /@ dataGrouped;
        sol = Map[Interpolation[#, InterpolationOrder -> 1]&, dataGrouped];
        res = 1. / diameter * Quiet @ NIntegrate[#[x], {x, 0., diameter}]& /@ sol;
        res
    ];


calROz[onedata_, yCoord_, diameter_, resREs_] :=
    Module[{sol, calData, dataGrouped, pos, res},
        calData = {#1, (1. - #8) * #5}& @@@ onedata;
        pos = findPosition[onedata, #]& /@ yCoord;
        dataGrouped = Pick[calData, #, 1]& /@ pos;
        dataGrouped = DeleteDuplicatesBy[#, First]& /@ dataGrouped;
        sol = Map[Interpolation[#, InterpolationOrder -> 1]&, dataGrouped];
        res = 1. / (1. - resREs[[#]]) / diameter * 
              Quiet @ NIntegrate[sol[[#]][x], {x, 0., diameter}]& /@ Range[Length@resREs];
        res
    ];


End[];


EndPackage[];


(* ::Example of package using:: *)
(*
The function "importNcal" read transient data of one operating 
condition, and calculate the cross-sectional averaged profiles. ;
	Input: address of the transient datas of one operating conditon. ;
	Output: an association of time averaged solids holdup, particle 
            velocity, and ozone concentration with y-coordinates. ;

importNcal[addr_] :=
 Block[{dir, rawdata, resREs, resRVs, resROz},
  SetDirectory[NotebookDirectory[]~StringJoin~addr];
  dir = Import["./"];
  Print[Length@dir];
  rawdata = Import[#] & /@ Select[dir, StringMatchQ[#, "*csv"] &];
  Print[TableForm[rawdata[[1, 1]], 
    TableHeadings -> {Automatic, None}]];
  
  rawdata = rawdata[[All, 2 ;; -1]];
  resREs = ParallelMap[calREs[#, yCoord, diameter] &, rawdata];
  resRVs = ParallelMap[calRVs[#, yCoord, diameter] &, rawdata];
  resROz = 
   ParallelMap[calROz[rawdata[[#]], yCoord, diameter, resREs[[#]]] &, 
    Range[Length@resREs]];
  
  Print["Calculation of " <> addr <> " complete."];
  
  Return[<|
    "es" -> Transpose@{Mean@resREs, yCoord},
    "vs" -> Transpose@{-1780. Mean@resRVs, yCoord},
    "oz" -> Transpose@{Mean@resROz, yCoord}
    |>];
  ]

  *)