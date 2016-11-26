function [ deltaV, isFound] = differential_corrector( init_state, t, yt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tolerance = 1e-6;
leftborder = 0-tolerance;
rightborder = 0+tolerance;

    if (yt(4) > leftborder && yt(4) < rightborder)  && (yt(6) > leftborder && yt(6) < rightborder)
        %deltaV = init_state - yt - yt/force_model(t,yt);
        isFound = true;
    else
        deltaV = (yt - (yt/force_model(t,yt))) - init_state;
        isFound = false;
    end

end

