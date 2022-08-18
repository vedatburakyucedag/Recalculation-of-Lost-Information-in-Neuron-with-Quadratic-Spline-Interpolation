function [v, w] = function_fhn_1n(a,b,c,I,h,N,t,v,w,aa,f,wet)
   
% define function handles
fv=@(t,v,w)  wet*(c * (v - w + I - (v^3) / 3) + (aa/(2*pi*f))*cos(2*pi*f*t));
fw=@(t,v,w)  wet*((v - b*w + a)/c);

    for i=1:N

        t(i+1) = t(i) + h;

        k1v = fv(t(i)     ,v(i)          ,w(i)          );
        k1w = fw(t(i)     ,v(i)          ,w(i)          );

        k2v = fv(t(i)+h/2 , v(i)+h/2*k1v , w(i)+h/2*k1w );
        k2w = fw(t(i)+h/2 , v(i)+h/2*k1v , w(i)+h/2*k1w );

        k3v = fv(t(i)+h/2 , v(i)+h/2*k2v , w(i)+h/2*k2w );
        k3w = fw(t(i)+h/2 , v(i)+h/2*k2v , w(i)+h/2*k2w );

        k4v = fv(t(i)+h   , v(i)+h  *k3v , w(i)+h  *k3w );
        k4w = fw(t(i)+h   , v(i)+h  *k3v , w(i)+h  *k3w );

        v(i+1) = v(i) + h/6 * (k1v + 2*k2v + 2*k3v + k4v);
        w(i+1) = w(i) + h/6 * (k1w + 2*k2w + 2*k3w + k4w);

    end
    
end

