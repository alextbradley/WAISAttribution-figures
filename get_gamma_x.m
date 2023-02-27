function gamma_x = get_gamma_x(x, slr, gamma)
%return the intersections of the curve slr = slr(gamma) with the curve slr = x
    [x0,~] = intersections(gamma, slr,[min(gamma), max(gamma)],x*[1,1]);
    gamma_x = x0;
end