(* ::Package:: *)

BeginPackage["exp`"];


calExpREs::usage = "Function to calculate the cross-sectional average solids holdup";
calExpRVs::usage = "Function to calculate the cross-sectional average solids velocity";
calExpROz::usage = "Function to calculate the cross-sectional average ozone concentration";


Begin["`Private`"];


calExpREs[data_, expPosRadial_, radium_] :=
    Block[{s},
        s = Interpolation[Transpose[{expPosRadial * radium, radium * 
            expPosRadial * data}], InterpolationOrder -> 1];
        2. / radium^2 * NIntegrate[s[x], {x, 0., radium}]
    ]


calExpRVs[i_, expDataUs_, expDataEs_, expMeanEs_, expPosRadial_, radium_] :=
    Block[{s},
        s = Interpolation[Transpose[{expPosRadial * radium, radium * 
            expPosRadial * (expDataUs * expDataEs)[[i]]}], InterpolationOrder -> 1];
        2. / radium^2 * NIntegrate[s[x], {x, 0., radium}]
    ]


calExpROz[i_, expDataOz_, expDataEs_, expMeanEs_, expPosRadial_, radium_] :=
    Block[{s},
        s = Interpolation[Transpose[{expPosRadial * radium, radium * 
            expPosRadial * (expDataOz * (1. - expDataEs))[[i]]}], InterpolationOrder -> 1];
        1. / (1. - expMeanEs[[i]]) * 2. / radium^2 * NIntegrate[s[x], {x, 0., radium}]
    ]


End[];


EndPackage[];


(* ::Example of package using:: *)
(*
expMeanEs5300 = 
  Quiet@Map[calExpREs[#, PositionRadial, radium] &, expData5300["es"]];

expMeanVs5300 = 
  Quiet@calExpRVs[#, expData5300["vs"], expData5300["es"], 
      expMeanEs5300, PositionRadial, radium] & /@ 
   Range[Length@PositionRiserWang];
   
expMeanOz5300 = 
  Quiet@calExpROz[#, expData5300["oz"], expData5300["es"], 
      expMeanEs5300, PositionRadial, radium] & /@ 
   Range[Length@PositionRiserWang];
*)