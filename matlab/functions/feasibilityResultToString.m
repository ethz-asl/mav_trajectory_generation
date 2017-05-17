function [ result_string ] = feasibilityResultToString( result_int )
%FEASIBILITYRESULTTOSTRING According to feasibility_base.cpp.
switch result_int
case 0
  result_string = 'Feasible';
case 1
  result_string = 'Indeterminable';
case 2
  result_string = 'ThrustHigh';
case 3
  result_string = 'ThrustLow';
case 4
  result_string = 'VelocityHigh';
case 5
  result_string = 'RollPitchRatesHigh';
case 6
  result_string = 'YawRatesHigh';
case 7
  result_string = 'YawAccHigh';
otherwise
  result_string = 'Unknown!';
end

