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
(* initial ozone mass fraction in Cheng-xiu's experiments, 100ppmv *)
Y0Li = 2.1649*^-5;
(* reaction constant in Cheng-xiu's experiments *)
krWang = 49.21 / av;
krLi = 4. / av;

(* 
------------------------------------------------------------------ 
-                           Functions 
------------------------------------------------------------------ 
*)

(* Drag Model *)

calKGdsp =
    Compile[
        {{epss, _Real}, {Uslip, _Real}}
        ,
        Module[{Cd, epsg, numRes, ksg, UslipC},
            If[epss == 0. || Uslip == 0.,
                Return[0.]
            ];
            UslipC = Abs @ Uslip;
            epsg = 1. - epss;
            numRes = epsg * rhog * UslipC * dp / mug;
            Cd = 24. / numRes (1. + .15 numRes^0.687);
            epsg = 1. - epss;
            If[epsg > 0.8,
                (* Wen-Yu *)
                ksg = .75 * Cd * (epss * epsg * rhog * UslipC) / dp * epsg ^ -2.65
                ,
                (* Ergun *)
                ksg = 150. * (mug * epss^2.) / (epsg * dp^2.) + 1.75 * (rhog * epss * UslipC) / dp
            ];
            If[Uslip > 0.,
                ksg
                ,
                -ksg
            ]
        ]
    ];

(* Sherwood Number *)

calKm =
    Compile[
        {{Uslip, _Real}, {epss, _Real}, {Xab, _Real}}
        ,
        Module[{Shlocal, Rep, Sc, epsg, tau, b, c, d, Sheff, Ug, kd},
            
            If[epss == 0. || Uslip == 0.,
                Return[0.]
            ];
            epsg = 1. - epss;
            tau = 1.2;
            b = .7;
            c = .5;
            d = 1. / 3.;
            (**)
            Ug = Abs[Uslip] * epsg;
            Rep = (rhog * Ug * dp) / mug;
            Sc = mug / (rhog * Dozone);
            Shlocal = (2. * epsg / tau) / (1. - CubeRoot[(1 - epsg) *
                 Xab]) + b * (Rep / epsg) ^ c * Sc^d;
            Sheff = (Rep^2 * Sc^2) / (12. * (1. - epsg) * Xab * epsg 
                / tau) * (Sqrt[1. + (24. * (1. - epsg) * Xab * epsg / tau) / (Rep^2 *
                 Sc^2) * Shlocal] - 1.);
            kd = Sheff * Dozone / dp;
            kd
        ]
    ];

(* Loaction and filter box function *)

(* find point's positions inside one frame *)

posInFrame2D[dataCal_, {{yLow_, yUp_}, {xLow_, xUp_}}] :=
    Boole[yLow <= #3 <= yUp && xLow <= #2 <= xUp]& @@@ dataCal;

(* find all filter box ranges, n is n times of dp *)

(* the output is {{yLow_,yUp_},{xLow_,xUp_}} *)
funcBox2D[n_, dataCal_] :=
    Module[{ylb, yub, xlb, xub, binWidth, ii, output},
        binWidth = n * dp;
        ylb = Range[0, Max[#3& @@@ dataCal] - binWidth, .5 binWidth];
        yub = ylb + binWidth;
        xlb = Range[0, Max[#2& @@@ dataCal] - binWidth, .5 binWidth];
        xub = xlb + binWidth;
        output = {};
        For[ii = 1, ii <= Length @ ylb, ii++,
            output = Append[output, {{ylb[[ii]], yub[[ii]]}, {xlb[[#]],
                 xub[[#]]}}& /@ Range[1, Length @ xlb]]
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

cal1FilterBox[dataCal_, bin_, nt_, daX_, daY_, daEs_, daP_, daUgxEg_,
     daUsxEs_, daUgyEg_, daUgyEg_, daUsyEs_, daOzEg_, daUslipDrag_, daOzEs_,
     dakm_] :=
    Block[{filterBoxPos, avX, avY, avEs, avP, avUgX, cUgX, avUsX, cUsX,
         avUslip, cUgY, avUgY, cUsY, avUsY, avUslipY, cOz, avOz, cOzEs, ckm, 
        avkm, cBetafy, avBetafy, avBetay, coefHdy, coefHr, coefHm},
        filterBoxPos = posInFrame2D[dataCal, bin];
        (************************************)
        avX = Mean[Pick[daX, filterBoxPos, 1]];
        avY = Mean[Pick[daY, filterBoxPos, 1]];
        avEs = Mean[Pick[daEs, filterBoxPos, 1]];
        avP = Mean[Pick[daP, filterBoxPos, 1]];
        (************************************)
        cUgX = Pick[daUgxEg, filterBoxPos, 1];
        avUgX = 1. / (1. - avEs) * Mean[cUgX];
        cUsX = Pick[daUsxEs, filterBoxPos, 1];
        avUsX = 1. / avEs * Mean[cUsX];
        cUgY = Pick[daUgyEg, filterBoxPos, 1];
        avUgY = 1. / (1. - avEs) * Mean[cUgY];
        cUsY = Pick[daUsyEs, filterBoxPos, 1];
        avUsY = 1. / avEs * Mean[cUsY];
        avUslip = Sqrt[(avUgX - avUsX) ^ 2 + (avUgY - avUsY) ^ 2];
        avUslipY = avUgY - avUsY;
        (************************************)
        cOz = Pick[daOzEg, filterBoxPos, 1];
        avOz = 1. / (1. - avEs) * Mean[cOz];
        (************************************)
        (************************************)
        cBetafy = Pick[daUslipDrag, filterBoxPos, 1];
        avBetafy = Mean[cBetafy] / avUslipY;
        avBetay = calKGdsp[avEs, avUslipY];
        coefHdy = avBetafy / avBetay;
        (************************************)
        cOzEs = Pick[daOzEs, filterBoxPos, 1];
        coefHr = Mean[cOzEs] / (avOz avEs);
        avkm = calKm[avUslip, avEs, 1.];
        ckm = Pick[dakm, filterBoxPos, 1];
        coefHm = Mean[cOzEs * ckm] / (avOz * avEs * avkm);
        (************************************)
        {nt, avX, avY, avEs, avP, avUgX, avUsX, avUgY, avUsY, avUslipY,
             avUslip, avOz, avkm, avBetay, coefHdy, coefHr, coefHm}
    ];

(* Calculate Results in One Filter Box Size & Export Table *)

cal1BoxSize[numData_, dataCal_, nt_, daX_, daY_, daEs_, daP_, daUgxEg_,
     daUsxEs_, daUgyEg_, daUgyEg_, daUsyEs_, daOzEg_, daUslipDrag_, daOzEs_,
     dakm_] :=
    Block[{boxBin, res},
        boxBin = funcBox2D[nt, dataCal];
        Print["Filter box bin calculation complete."];
        res = ParallelMap[cal1FilterBox[dataCal, #, nt, daX, daY, daEs,
             daP, daUgxEg, daUsxEs, daUgyEg, daUgyEg, daUsyEs, daOzEg, daUslipDrag,
             daOzEs, dakm]&, boxBin];
        Print["Results calculation complete."];
        (* Export Results *)
        Export["table\\coef-" ~ StringJoin ~ ToString[nt] ~ StringJoin ~ "dp-" ~ StringJoin ~ ToString[dir[[numData]]], res];
        Print["---------- Results Export Complete ----------"];
    ];

(* Calculate Parameters for One Data Set *)

cal1DataSet[numData_] :=
    Block[
        {dataCal, daX, daY, daP, daOz, daUgx, daUgy, daUsx, daUsy, daEs,
             daUslip, daUgxEg, daUsxEs, daUgyEg, daUsyEs, daOzEg, daOzEs, dakm, daDrag,
             daUslipDrag}
        ,
        (* Data Reading *)
        Print["---------------------------------------------"];
        Print["---------------- NEW DATA SET ---------------"];
        Print["-------------------- ", numData, "/", Length @ dir, " --------------------"
            ];
        Print["---------------------------------------------"];
        Print["The data set is " ~ StringJoin ~ ToString[dir[[numData
            ]]]];
        dataCal = Import["data\\" ~ StringJoin ~ ToString[dir[[numData
            ]]]][[2 ;; -1]];
        Print["Data reading complete."];
        (* Calculate parameter set 1 *)
        daX = #2& @@@ dataCal;
        daY = #3& @@@ dataCal;
        daOz = #4& @@@ dataCal;
        daP = #5& @@@ dataCal;
        daUgx = #6& @@@ dataCal;
        daUgy = #7& @@@ dataCal;
        daUsx = #8& @@@ dataCal;
        daUsy = #9& @@@ dataCal;
        daEs = #10& @@@ dataCal;
        Print["Parameter set 1 calculation complete."];
        (* Calculate parameter set 2 *)
        daUslip = Sqrt[(daUgx - daUsx) ^ 2 + (daUgy - daUsy) ^ 2];
        daUgxEg = daUgx * (1. - daEs);
        daUsxEs = daUsx * daEs;
        daUgyEg = daUgy * (1. - daEs);
        daUsyEs = daUsy * daEs;
        daOzEg = daOz * (1. - daEs);
        daOzEs = daOz * daEs;
        daDrag = ParallelMap[calKGdsp[daEs[[#]], daUgy[[#]] - daUsy[[
            #]]]&, Range[Length @ daEs]];
        daUslipDrag = (daUgy - daUsy) * daDrag;
        dakm = ParallelMap[calKm[daUslip[[#]], daEs[[#]], 1.]&, Range[
            Length @ daEs]];
        Print["Parameter set 2 calculation complete."];
        (* Calculate Coefficients *)
        cal1BoxSize[numData, dataCal, #, daX, daY, daEs, daP, daUgxEg,
             daUsxEs, daUgyEg, daUgyEg, daUsyEs, daOzEg, daUslipDrag, daOzEs, dakm
            ]& /@ range;
    ];

(* 
------------------------------------------------------------------ 
-                             Export 
------------------------------------------------------------------ 
*)

range = {27, 36, 48};

Map[cal1DataSet, Range[Length @ dir]];

CloseKernels[];

(* CMD line running code: wolframscript -file '.\cal_filter_2D.wl' *)