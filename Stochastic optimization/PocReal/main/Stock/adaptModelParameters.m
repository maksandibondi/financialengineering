function [inputStruct] = adaptModelParameters( inputStruct, parameterToAdapt)
%ADAPTMODELPARAMETERS Summary of this function goes here
%   Detailed explanation goes here

if strcmp(parameterToAdapt,'epsilon')
    inputStruct.epsilon = inputStruct.epsilon + 0.02;   
elseif strcmp(parameterToAdapt,'weight')
    inputStruct.concentration_weights = inputStruct.concentration_weights*0.75;   
end

