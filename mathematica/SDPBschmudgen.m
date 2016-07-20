(* ::Package:: *)

(* Load functions from SDPB.m *)
<<"SDPB.m"

(* Sample points inside the interval should be ignored *)
rescaledLaguerreSamplePoints[n_, interval_] := Module[{npts=0, k=0, samplePoints={}, pt},
	If[Length[interval]==0, Return[rescaledLaguerreSamplePoints[n]];];
	While[npts < n,
	pt = SetPrecision[\[Pi]^2 (-1+4k)^2/(-64n Log[rhoCrossing]), prec];
    If[!(interval[[1]]< SetPrecision[\[Pi]^2 (-1+4k)^2/(-64n Log[rhoCrossing]), prec]<interval[[2]]),
		AppendTo[samplePoints,pt];npts++;]; k++;];
	Return[samplePoints];	
];

(* Write schmudgenInterval XML node when PositiveMatrixWithPrefactor has interval specified *)
ClearAll[WriteBootstrapSDP];
WriteBootstrapSDP[file_, SDP[objective_, normalization_, positiveMatricesWithPrefactors_]] := Module[
    {
        stream = OpenWrite[file],
        node, real, int, vector, polynomial,
        polynomialVector, polynomialVectorMatrix,
        affineObjective, polynomialVectorMatrices
    },

    (* write a single XML node to file.  children is a routine that writes child nodes when run. *)
    node[name_, children_] := (
        WriteString[stream, "<", name, ">"];
        children[];
        WriteString[stream, "</", name, ">\n"];
    );

    real[r_][] := WriteString[stream, nf[r]];
    int[i_][] := WriteString[stream, i];
    vector[v_][] := Do[node["elt", real[c]], {c, v}];
    polynomial[p_][] := Do[node["coeff", real[c]], {c, safeCoefficientList[p,x]}];
    polynomialVector[v_][] := Do[node["polynomial", polynomial[p]], {p, v}];

    polynomialVectorMatrix[PositiveMatrixWithPrefactor[prefactor_, m_, interval_:{}]][] := Module[
        {degree = Max[Exponent[m, x]], samplePoints, sampleScalings, bilinearBasis},
		samplePoints   = rescaledLaguerreSamplePoints[degree + 1, interval];
        sampleScalings = Table[prefactor /. x -> a, {a, samplePoints}];
        bilinearBasis  = orthogonalPolynomials[prefactor, Floor[degree/2]];
        node["rows", int[Length[m]]];
        node["cols", int[Length[First[m]]]];
        node["elements", Function[
            {},
            Do[node[
                "polynomialVector",
                polynomialVector[reshuffleWithNormalization[normalization,pv]]],
               {row, m}, {pv, row}]]];
        node["samplePoints", vector[samplePoints]];
        node["sampleScalings", vector[sampleScalings]];
        node["bilinearBasis", polynomialVector[bilinearBasis]];
		If[Length[interval]>0, node["schmudgenInterval", vector[interval]];]

    ];

    node["sdp", Function[
        {},
        node["objective", vector[reshuffleWithNormalization[normalization, objective]]];
        node["polynomialVectorMatrices", Function[
            {},
            Do[node["polynomialVectorMatrix", polynomialVectorMatrix[pvm]], {pvm, positiveMatricesWithPrefactors}];
        ]];
    ]];                                          

    Close[stream];
];
