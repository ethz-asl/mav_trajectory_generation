function [ result_string ] = feasibilityResultToString( result_int )
%FEASIBILITYRESULTTOSTRING According to feasibility_base.cpp.
switch result_int
case 0
  result_string = "Feasible";
case 1
  result_string = "Indeterminable";
case 2
  result_string = "InfeasibleThrustHigh";
case 3
  result_string = "InfeasibleThrustLow";
case 4
  result_string = "InfeasibleVelocity";
case 5
  result_string = "InfeasibleRollPitchRates";
case 6
  result_string = "InfeasibleYawRates";
case 7
  result_string = "InfeasibleYawAcc";
otherwise
  result_string = "Unknown!";
end

