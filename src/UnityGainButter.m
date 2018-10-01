(* ::Package:: *)

(* :Title: Unity-Gain Butterworth Filters*)
(* :Author: David E. McKnight*)
(* :Version: May 2018*)
(* :Summary:
    A package for generating arbitrarily-ordered Unity-Gain Butterworth filters with specified cutoff frequency
 *)
 
BeginPackage["UnityGainButterworth`"]

(**Low-Pass filter generator declaration**)
LPUnityButter::usage = "LPUnityButter[n, wc, var] creates an nth-order, low-pass, unity-gain Butterworth filter with cutoff frequency wc (in radians) under the given variable var";

Begin["Private`"]
(**Low-Pass filter generator definition**)
LPUnityButter[n_, wc_, var_]:=
Module[
{n0=n, wc0=wc, var0=var},
(**Find the required denominator**)
denom=getButterworthDenom[n0, wc0, var];
(**Find the required numerator**)
num = Function[k, denom/.var->k][0];

(**Simplify expression**)
expression = num/denom;
newNum = Numerator[expression];
newDenom = Denominator[expression];

(**Remove imaginary coefficients**)
denomList = Reverse[CoefficientList[newDenom, var]];
numList = Reverse[CoefficientList[newNum,var]];
fixedDenomList=Re/@denomList;
fixedNumList=Re/@numList;
fixedDenom2=Expand[FromDigits[fixedDenomList, var]];
fixedNum2=Expand[FromDigits[fixedNumList, var]];

(**Return the final expression**)
fixedNum2/fixedDenom2
]

(**Retrieve the denominator of an nth-order Butterworth filter with the given corner frequency**)
getButterworthDenom[n_, wc_, var_]:=Expand[factorForm[leftHalfPlaneRoots[1/magSimpButterSquared[n,wc, var], var],var]] 
(**Returns the magnitude of a primitive nth-order Butterworth filter with the given cutoff frequency**)
magSimpButterSquared[n_, wc_, var_]:=1/(1+(-1)^n*var^(2n)/wc^(2n));
(**Puts roots into the form (s-root1)(s-root2)...(s-lastRoot)**)
factorForm[rootsList_, var_]:=Expand[Times@@(var-#&/@rootsList)];
(**Grabs the left-half-plane roots**)
leftHalfPlaneRoots[rootsList_, var_]:=Select[rootsToList[rootsList,var],Re[#]<=0&]
(**Creates a list of the roots**)
rootsToList[a_,var_]:=rootsToDisjunction[a,var][[All,-1]];
rootsToDisjunction[a_, var_]:=Apply[List,ComplexExpand[Roots[a==0,var]]];
End[]

EndPackage[]
