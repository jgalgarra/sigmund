function [y] = pruebaprob (pob,vp)
    a = pob;
    rper = rp(vp,365)
    for i=1:365
        a = a + sign(rper)*binornd(a,abs(rper))
    end
    y = a;