LaunchKernels[];
Print["The total number of cores is ", Length @ Kernels[]];
dir = Import["data"];
Print["The total number of data is ", Length @ dir];

(* 
------------------------------------------------------------------ 
-                        Basic Informations 
------------------------------------------------------------------ 
*)

rhog = 1.225;
dp = 70.*^-6;
mug = 1.7894*^-5;
rhos = 1780.;
av = 6. / dp;
Dozone = 1.48535*^-5;
(* initial ozone mass fraction in Cheng-xiu's experiments, 100ppmv *)
Y0Wang = 1.6653*^-4;
(* initial ozone mass fraction in Dong-bin's experiments, 13ppmv *)
Y0Li = 2.1649*^-5;
(* reaction constant in Cheng-xiu's experiments *)
krWang = 49.2 / av;
krLi = 4. / av;

(* 
------------------------------------------------------------------ 
-                           Functions 
------------------------------------------------------------------ 
*)

(* Drag Model *)

calKGdsp = Compile[
    {{epss, _Real}, {Uslip, _Real}}
    ,
    Module[{Cd, epsg, numRes, ksg, UslipC},
        If[epss == 0. || Uslip == 0.,
            Return[0.]
        ];
        UslipC = Abs @ Uslip;
        epsg = 1. - epss;
        numRes = epsg * rhog * UslipC * dp / mug;
        Cd = 24. / numRes (1. + .15*numRes^0.687);
        epsg = 1. - epss;
        If[epsg > 0.8,
            (* Wen-Yu *)
            ksg = .75 * Cd * (epss * epsg * rhog * UslipC) / dp * epsg^-2.65
            ,
            (* Ergun *)
            ksg = 150. * (mug * epss^2.) / (epsg * dp^2.) + 1.75 * (rhog * epss * UslipC) / dp
        ];
        If[Uslip > 0., ksg , -ksg]
    ]
];

(* Loaction and filter box function *)

(* find point's positions inside one frame *)
posInFrame[dataCoord_, {{yLow_, yUp_}, {xLow_, xUp_}, {zLow_, zUp_}}] :=
    Block[{x, y, z},
        x = dataCoord[[All, 1]];
        y = dataCoord[[All, 2]];
        z = dataCoord[[All, 3]];
        Boole[yLow < y[[#]] < yUp && xLow < x[[#]] < xUp && zLow < z[[#]] < zUp] & /@ Range[Length@x]
    ];
(* find all filter box ranges, n is n times of dp *)

(* the output is {{yLow_,yUp_},{xLow_,xUp_},{zLow_,zUp_}} *)
funcBox[n_, dataCal_] := 
    Module[{ylb, yub, xlb, xub, zlb, zub, binWidth, binError, ii, jj, output},
        binWidth = n*dp;
        binError = 3.3 dp/3;
        ylb = Range[0, Max[#2 & @@@ dataCal] - binWidth, 0.5*binWidth]-binError;
        yub = ylb + binWidth+binError;
        xlb = Range[0, Max[#1 & @@@ dataCal] - binWidth, 0.5*binWidth]-binError;
        xub = xlb + binWidth+binError;
        zlb = Range[0, Max[#3 & @@@ dataCal] - binWidth, 0.5*binWidth]-binError;
        zub = zlb + binWidth+binError;

        output = {};
        For[jj = 1, jj <= Length@ylb, jj++,
            For[ii = 1, ii <= Length@zlb, ii++,
                output = 
                Append[output, {{ylb[[jj]], yub[[jj]]}, {xlb[[#]], 
                    xub[[#]]}, {zlb[[ii]], zub[[ii]]}} & /@ Range[1, Length@xlb]]
            ]
        ];
        Print["The filter box size is ", n, " times of dp, which is ",
            binWidth 100., " cm."];
        Print["The number of box is ", Length @ Flatten[output, 1]];
        Flatten[output, 1]
    ];

(* 
------------------------------------------------------------------
-                     Calculation Functions 
------------------------------------------------------------------ 
*)

(* Calculatee Parameters in One Box *)

cal1FilterBox[dataCoord_, bin_, nt_,
                daX_, daY_, daZ_,
                daEs_, daP_,
                daUgxEg_, daUsxEs_, 
                daUgyEg_, daUsyEs_, 
                daUgzEg_, daUszEs_, 
                daOzEg_, daUslipDrag_, daOzEs_] :=
    Block[{filterBoxPos,
            avX, avY, avZ,
            avEs, avP,
            avUgX, cUgX,
            avUsX, cUsX,
            avUgY, cUgY,
            avUsY, cUsY,
            avUgZ, cUgZ,
            avUsZ, cUsZ,
            avUslip, avUslipY,
            cOz, avOz, cOzEs,
            cBetafy, avBetafy, avBetay,
            coefHdy, coefHr},
        filterBoxPos = posInFrame[dataCoord, bin];
        (************************************)
        avX = Mean[Pick[daX, filterBoxPos, 1]];
        avY = Mean[Pick[daY, filterBoxPos, 1]];
        avZ = Mean[Pick[daZ, filterBoxPos, 1]];
        avEs = Mean[Pick[daEs, filterBoxPos, 1]];
        avP = Mean[Pick[daP, filterBoxPos, 1]];
        (************************************)
        cUgX = Pick[daUgxEg, filterBoxPos, 1];
        avUgX = Mean[cUgX]/(1. - avEs);
        (**)
        cUsX = Pick[daUsxEs, filterBoxPos, 1];
        avUsX = Mean[cUsX]/avEs;
        (**)
        cUgY = Pick[daUgyEg, filterBoxPos, 1];
        avUgY = Mean[cUgY]/(1. - avEs);
        (**)
        cUsY = Pick[daUsyEs, filterBoxPos, 1];
        avUsY = Mean[cUsY]/avEs;
        (**)
        cUgZ = Pick[daUgzEg, filterBoxPos, 1];
        avUgZ = Mean[cUgZ]/(1. - avEs);
        (**)
        cUsZ = Pick[daUszEs, filterBoxPos, 1];
        avUsZ = Mean[cUsZ]/avEs;
        (**)
        avUslip = Sqrt[(avUgX - avUsX)^2 + (avUgY - avUsY)^2 + (avUgZ - avUsZ)^2];
        avUslipY = avUgY - avUsY;
        (************************************)
        cOz = Pick[daOzEg, filterBoxPos, 1];
        avOz = Mean[cOz]/(1. - avEs);
        (************************************)
        (************************************)
        cBetafy = Pick[daUslipDrag, filterBoxPos, 1];
        avBetafy = Mean[cBetafy] / avUslipY;
        avBetay = calKGdsp[avEs, avUslipY];
        coefHdy = avBetafy / avBetay;
        (************************************)
        cOzEs = Pick[daOzEs, filterBoxPos, 1];
        coefHr = Mean[cOzEs] / (avOz * avEs);
        (************************************)
        {nt, 
         avX, avY, avZ, 
         avEs, avP, 
         avUgX, avUsX, avUgY, avUsY, avUgZ, avUsZ, 
         avUslipY, avUslip, 
         avOz, 
         avBetay,
         coefHdy, coefHr}
    ];

(* Calculate Results in One Filter Box Size & Export Table *)

cal1BoxSize[numData_, dataCoord_, nt_, 
            daX_, daY_, daZ_, daEs_, daP_, 
            daUgxEg_, daUsxEs_, 
            daUgyEg_, daUsyEs_, 
            daUgzEg_, daUszEs_, 
            daOzEg_, daUslipDrag_, daOzEs_] :=
    Block[{boxBin, res},
        boxBin = funcBox[nt, dataCoord];
        Print["The number of cell in one box is ", Total[posInFrame[dataCoord, boxBin[[1]]]]];
        Print["Filter box bin calculation complete."];
        res = ParallelMap[cal1FilterBox[dataCoord, #, nt, 
                                        daX, daY, daZ, 
                                        daEs, daP, 
                                        daUgxEg, daUsxEs, 
                                        daUgyEg, daUsyEs, 
                                        daUgzEg, daUszEs, 
                                        daOzEg, daUslipDrag, daOzEs]&, boxBin];
        Print["Results calculation complete."];
        (* Export Results *)
        Export["table\\coef-" ~StringJoin~ ToString[nt] ~StringJoin~ "dp-" ~StringJoin~ ToString[dir[[numData]]], res];
        Print["---------- Results Export Complete ----------"];
    ];

(* Calculate Parameters for One Data Set *)

cal1DataSet[numData_] :=
    Block[{dataCal,
            daX, daY, daZ, dataCoord,
            daP, daOz,
            daUgx, daUgy, daUgz,
            daUsx, daUsy, daUsz,
            daEs, daUslip,
            daUgxEg, daUsxEs,
            daUgyEg, daUsyEs,
            daUgzEg, daUszEs,
            daOzEg, daOzEs, daDrag, daUslipDrag}
        ,
        (* Data Reading *)
        Print["---------------------------------------------"];
        Print["---------------- NEW DATA SET ---------------"];
        Print["-------------------- ", numData, "/", Length @ dir, " --------------------"];
        Print["---------------------------------------------"];
        Print["The data set is " ~StringJoin~ ToString[dir[[numData]]]];
        dataCal = Import["data\\" ~StringJoin~ ToString[dir[[numData]]]][[2 ;; -1]];
        Print["Data reading complete."];
        (* Calculate parameter set 1 *)
        daX = dataCal[[All, 2]];
        daY = dataCal[[All, 3]];
        daZ = dataCal[[All, 4]];
        dataCoord = Transpose@{daX, daY, daZ};
        daOz = dataCal[[All, 5]];
        daP = dataCal[[All, 6]];
        daUgx = dataCal[[All, 7]];
        daUgy = dataCal[[All, 8]];
        daUgz = dataCal[[All, 9]];
        daUsx = dataCal[[All, 10]];
        daUsy = dataCal[[All, 11]];
        daUsz = dataCal[[All, 12]];
        daEs = dataCal[[All, 13]];
        Print["Parameter set 1 calculation complete."];
        (* Calculate parameter set 2 *)
        daUslip = Sqrt[(daUgx - daUsx)^2 + (daUgy - daUsy)^2 + (daUgz - daUsz)^2];
        daUgxEg = daUgx*(1. - daEs);
        daUsxEs = daUsx*daEs;
        daUgyEg = daUgy*(1. - daEs);
        daUsyEs = daUsy*daEs;
        daUgzEg = daUgz*(1. - daEs);
        daUszEs = daUsz*daEs;
        daOzEg = daOz*(1. - daEs);
        daOzEs = daOz*daEs;
        daDrag = ParallelMap[calKGdsp[daEs[[#]], daUgy[[#]] - daUsy[[#]]] &, Range[Length@daEs]];
        daUslipDrag = (daUgy - daUsy)*daDrag;
        Print["Parameter set 2 calculation complete."];
        (* Calculate Coefficients *)
        cal1BoxSize[numData, dataCoord, #, 
                    daX, daY, daZ, 
                    daEs, daP, 
                    daUgxEg, daUsxEs, 
                    daUgyEg, daUsyEs, 
                    daUgzEg, daUszEs, 
                    daOzEg, daUslipDrag, daOzEs
                ]& /@ range;
    ];

(* 
------------------------------------------------------------------ 
-                             Export 
------------------------------------------------------------------ 
*)

range = {30, 40, 50, 80, 100};

Map[cal1DataSet, Range[Length @ dir]]//AbsoluteTiming

CloseKernels[];

(* CMD line running code: wolframscript -file '.\cal_filter_3D.wl' *)