(*
    This codes calculates the gradients and generates table contain the
    gradients and other parameter in the raw dataset.
    The input datasets are generated from codes: 'cal_filter.wl'.

    The input dataset's structure should be:
    {nt
    , avX, avY
    , avEs, avP
    , avUgX, avUsX, avUgY, avUsY
    , avUslipY, avUslip
    , avOz
    , avkm, avBetay
    , coefHdy, coefHr, coefHm}
*)
(* 
    -----------------------------------------------------------------
    -                         Dataset Import
    -----------------------------------------------------------------
*)
dir = Import["table"];
Print["The total number of data is ", Length @ dir];
dataList = Import[dir];
listPos = Boole[StringSplit[#, "."][[-1]] == "csv"]& /@ dataList
dataList = Pick[dataList, listPos, 1];
dataTrain = Import[ToString[dir] ~ StringJoin ~ "/" ~ StringJoin ~ 
    ToString[#]]& /@ dataList;
Print[Total[Differences[Length[#]& /@ dataTrain[[All, 1]]]], 
    " = 0: the structures of datasets are the same. "];

(* 
    ------------------------------------------------------------------
    -                            Functions
    ------------------------------------------------------------------ 
*)
(* Evaluate derivate of a interpolation function *)
calGrad2D[itpltF_, x_] := Evaluate[D[itpltF[x], x]];

(* Tidy up the coordinates of one time dataset *)
funcReplaceCoordinate[data1Time_, precision_] :=
    Block[{coordinate, dataReplace},
        coordinate = {SetPrecision[#2& @@@ data1Time, precision]} 
            ~ Join ~ {SetPrecision[#3& @@@ data1Time, precision]};
        dataReplace = ReplacePart[data1Time[[#]], 2 -> coordinate[[1,
             #]]]& /@ Range[Length @ data1Time];
        dataReplace = ReplacePart[dataReplace[[#]], 3 -> coordinate[[
            2, #]]]& /@ Range[Length @ data1Time]
    ];

(* 
    Check the coordinates after tidy up 
*)
checkCoor[dataPar_] :=
    Block[{},
        Total[Total[(#2& @@@ dataPar[[#[[1]]]]) - (#2& @@@ dataPar[[#[[2]]]])
            ]& /@ Table[{RandomInteger[{1, Length @ dataPar}], 
            RandomInteger[{1, Length @ dataPar}]}, {i, 1, 10000}]] 
        + Total[Total[Differences[#3& @@@ #]]& /@ dataPar]
    ];

(* 
    input partition dataset at one time and give a table of gradient 
    for the dataset 
*)
gradPar[dataPar_, width_] := Block[{radialGrad, axialGrad, dataGrad},
    Print[
        "--------------------------------------------------
        The gradients include: Es, P, UslipY, Uslip, Oz.
        The gradient directions are: radial and axial. 
        --------------------------------------------------"];
    dataGrad = {};
    (*-----------------------------------------------------------*)
    (* Gradient of solids holdup #4 *)
    radialGrad = Block[{itplt, grad1},
        itplt = Interpolation[{#2, #4} & @@@ #, Method -> "Spline", 
                InterpolationOrder -> 3];
        grad1 = calGrad2D[itplt, #2] & @@@ #
    ] & /@ dataPar;
    AppendTo[dataGrad, Flatten@radialGrad];
   
    axialGrad = Block[{itplt},
        itplt = Interpolation[{#3, #4} & @@@ dataPar[[All, #]]];
        calGrad2D[itplt, #] & /@ (#3 & @@@ dataPar[[All, #]])
    ] & /@ Range[width];
    AppendTo[dataGrad, Flatten@Transpose[axialGrad]];
    (*-----------------------------------------------------------*)
    (* Gradient of Pressure #5 *)
    radialGrad = Block[{itplt, grad1},
        itplt = Interpolation[{#2, #5} & @@@ #, Method -> "Spline", 
                InterpolationOrder -> 3];
        grad1 = calGrad2D[itplt, #2] & @@@ #
    ] & /@ dataPar;
    AppendTo[dataGrad, Flatten@radialGrad];
   
    axialGrad = Block[{itplt},
        itplt = Interpolation[{#3, #5} & @@@ dataPar[[All, #]]];
        calGrad2D[itplt, #] & /@ (#3 & @@@ dataPar[[All, #]])
    ] & /@ Range[width];
    AppendTo[dataGrad, Flatten@Transpose[axialGrad]];
    (*-----------------------------------------------------------*)
    (* Gradient of slip velocity Y #10 *)
    radialGrad = Block[{itplt, grad1},
        itplt = Interpolation[{#2, #10} & @@@ #, Method -> "Spline", 
                InterpolationOrder -> 3];
        grad1 = calGrad2D[itplt, #2] & @@@ #
    ] & /@ dataPar;
    AppendTo[dataGrad, Flatten@radialGrad];
   
    axialGrad = Block[{itplt},
        itplt = Interpolation[{#3, #10} & @@@ dataPar[[All, #]]];
        calGrad2D[itplt, #] & /@ (#3 & @@@ dataPar[[All, #]])
    ] & /@ Range[width];
    AppendTo[dataGrad, Flatten@Transpose[axialGrad]];
    (*-----------------------------------------------------------*)
    (* Gradient of slip velocity #11 *)
    radialGrad = Block[{itplt, grad1},
        itplt = Interpolation[{#2, #11} & @@@ #, Method -> "Spline", 
                InterpolationOrder -> 3];
        grad1 = calGrad2D[itplt, #2] & @@@ #
    ] & /@ dataPar;
    AppendTo[dataGrad, Flatten@radialGrad];
   
    axialGrad = Block[{itplt},
       itplt = Interpolation[{#3, #11} & @@@ dataPar[[All, #]]];
       calGrad2D[itplt, #] & /@ (#3 & @@@ dataPar[[All, #]])
       ] & /@ Range[width];
    AppendTo[dataGrad, Flatten@Transpose[axialGrad]];
    (*-----------------------------------------------------------*)
    (* Gradient of ozone mass fraction #12 *)
    radialGrad = Block[{itplt, grad1},
        itplt = Interpolation[{#2, #12} & @@@ #, Method -> "Spline", 
                InterpolationOrder -> 3];
        grad1 = calGrad2D[itplt, #2] & @@@ #
    ] & /@ dataPar;
    AppendTo[dataGrad, Flatten@radialGrad];
   
    axialGrad = Block[{itplt},
        itplt = Interpolation[{#3, #12} & @@@ dataPar[[All, #]]];
        calGrad2D[itplt, #] & /@ (#3 & @@@ dataPar[[All, #]])
    ] & /@ Range[width];
    AppendTo[dataGrad, Flatten@Transpose[axialGrad]];
    (*------------------------- Result ---------------------------*)
  
    Print["The gradients generated successfully. "];
    dataGrad
];

(* 
    input dataset at one time, do coordinate tidy up, partition, 
    check, and give gradient table + tide-up dataset
*)
gradData[dataCal_] := Block[{width, dataPar, dataReplaced, table},
    Print["---------------------- START ---------------------"];
    Print["The filtering box size is ", dataCal[[1, 1]]];
    dataReplaced = funcReplaceCoordinate[dataCal, 8];
    width = Length @ DeleteDuplicates[#2& @@@ dataReplaced];
    Print["The radial width is ", width];
    dataPar = Partition[dataReplaced, width];
    Print["Check the coordinate tidy-up, '0' means good: ", 
        Boole[checkCoor[dataPar] != 0.]];
    table = gradPar[dataPar, width] ~ Join ~ Transpose[dataReplaced];
    Print["The data table generated successfully. "];
    Print["----------------------- END ----------------------"];
    table
];

(* 
    -----------------------------------------------------------------
    -                            Export
    -----------------------------------------------------------------
*)
Block[{},
    Print["---------------------- ", #, "/", 
        Length @ dataTrain, " ---------------------"];
    Export["tableGrad/grad-" 
        <> StringJoin[Riffle[StringSplit[dataList[[#]], "-"][[2 ;; -1]], "-"]], 
        gradData[dataTrain[[#]]]
    ] & /@ Range[Length @ dataTrain]
]