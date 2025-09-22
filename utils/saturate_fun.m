function saturation =  saturate_fun(x, x_max, x_min) 

    saturation = min(x_max, max(x_min,x)) ; % Saturation Function

end