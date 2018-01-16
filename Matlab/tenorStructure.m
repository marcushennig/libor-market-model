% Returns the tenor structure
function T = tenorStructure(tstart, tend, tenor)
    T = [];
    t = tstart;
    while t <= tend
        T = [ T  t ];
        t = t + calmonths(tenor); 
    end
    if(T(end) ~= tend)
        T = [ T tend ];
    end
end