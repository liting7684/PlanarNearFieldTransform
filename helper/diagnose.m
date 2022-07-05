function diagnose(fx,d,k0,kx,ky,Mc_padded,Nc_padded)
%DIAGNOSE 此处显示有关此函数的摘要
%   此处显示详细说明
    fx_new = fx.*exp(1i*sqrt(k0.^2-kx.^2-ky.^2)*d);
    E_aperture = ifft2(fx_new,Mc_padded,Nc_padded);
    surf(abs(E_aperture));
end

