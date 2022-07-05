%% 使用FFT的平面波谱与远场方向图计算
% variables:
% M,N: x,y方向天线单元的个数;M=43,N=31
% lambda: 电磁波波长
% Im: 半波振子波腹电流
% dx,dy: x,y方向天线单元的间距
% theta: theta角度范围精度
% phi: phi角度范围精度
% deltax,deltay: 在x轴和y轴上的采样间距
% d: 采样平面距离阵列平面的间距
% Mc,Nc: 计算机仿真x,y方向的采样点数目,取2的倍数
% 该方法为计算平面波谱推远场方向图
% Author:THY
%% 

function [data_nf2ff] = nf2ff_planar_fft(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,phi,d,Mc,Nc,padding)
    % ---------------------------------------------------------------------
    % 构建结果表
    % p是theta的长度，t是phi的长度
    p=length(theta);
    t=length(phi);
    data_nf2ff = table(zeros(p*t,1),zeros(p*t,1),zeros(p*t,3),zeros(p*t,1),zeros(p*t,1),zeros(p*t,1));
    data_nf2ff.Properties.VariableNames = {'theta','phi','E','Etheta','Ephi','Eabs'};
    % ---------------------------------------------------------------------
    % 计算波数k0
    k0=2*pi/lambda;
    % ---------------------------------------------------------------------
    Mc_padded=padding*Mc;    
    Nc_padded=padding*Nc;
    m=-Mc_padded/2:1:Mc_padded/2-1;
    n=-Nc_padded/2:1:Nc_padded/2-1;
    % ---------------------------------------------------------------------
    % kx=2*pi*m/(Mc*deltax)
    % ky=2*pi*n/(Nc*deltay)
    kx=2*pi*m/(Mc_padded*deltax);
    ky=2*pi*n/(Nc_padded*deltay);
    [ky_grid,kx_grid] = meshgrid(ky,kx);
    kz_grid = sqrt(k0^2-kx_grid.^2-ky_grid.^2);
    % ---------------------------------------------------------------------
    %球面远场波数矢量
    [theta_grid,phi_grid]=meshgrid(theta,phi);
    kx_grid_spherical = k0*sin(theta_grid).*cos(phi_grid);
    ky_grid_spherical = k0*sin(theta_grid).*sin(phi_grid);
    kz_grid_spherical = k0*cos(theta_grid);
    % ---------------------------------------------------------------------
    % 计算Ex
    Ex=0;
    for m=0:Mc-1
        for n=0:Nc-1
            XL1=(m-Mc/2)*deltax;
            YL2=(n-Nc/2)*deltay;
            Ex(m+1,n+1)=Ex_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im);%计算Ex的采样电场值
        end
    end
%     surf(abs(Ex));
%     surf(angle(Ex));
    % ---------------------------------------------------------------------
    % 计算Ey
    Ey=0;
    for m=0:Mc-1
        for n=0:Nc-1
            XL1=(m-Mc/2)*deltax;
            YL2=(n-Nc/2)*deltay;
            Ey(m+1,n+1)=Ey_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im,Mc,Nc,deltax,deltay);%计算Ey的采样电场值
        end
    end
%     surf(abs(Ey));
%     surf(angle(Ey));
    % ---------------------------------------------------------------------
    fx=ifftshift(ifft2(Ex,Mc_padded,Nc_padded));
    %mesh(kx/k0,ky/k0,(abs(fx)))
    fy=ifftshift(ifft2(Ey,Mc_padded,Nc_padded));
    %mesh(kx/k0,ky/k0,(abs(fy)))
    fz=-(fx.*kx_grid+fy.*ky_grid)./kz_grid;
    % ---------------------------------------------------------------------
    % 将fx平面波谱推到口面上
    % 参考文献:Digital Holographic Reconstruction and Filtering Method 
    % for Antenna Planar Near-Field Phase-Less Measurement
    diagnose(fx,d,k0,kx,ky,Mc_padded,Nc_padded);


    % ---------------------------------------------------------------------
    [kx,ky]=meshgrid(ky,kx);
    % Interpolate Modes in spherical coordinates
    fx_ff=interp2(kx,ky,abs(fx),kx_grid_spherical,ky_grid_spherical,'spline');
    fy_ff=interp2(kx,ky,abs(fy),kx_grid_spherical,ky_grid_spherical,'spline');
    fz_ff=interp2(kx,ky,abs(fz),kx_grid_spherical,ky_grid_spherical,'spline');
    % ---------------------------------------------------------------------
    % Far Field 
    r=100;
    C=1j*(k0*exp(-1j*k0*r))/(2*pi*r);
    Etheta=C*(fx_ff.*cos(phi_grid)+fy_ff.*sin(phi_grid));
    mesh(kx_grid_spherical,ky_grid_spherical,20*log(abs(Etheta)))
    Ephi=C*cos(theta_grid).*(-fx_ff.*sin(phi_grid)+fy_ff.*cos(phi_grid));
    mesh(kx_grid_spherical,ky_grid_spherical,20*log(abs(Ephi)))
    Ex = C*cos(theta_grid).*fx_ff;
    Ey = C*cos(theta_grid).*fy_ff;
    Ez = C*cos(theta_grid).*fz_ff;

    % ---------------------------------------------------------------------
    s = numel(theta_grid);
    data_nf2ff.Etheta = reshape(Etheta,s,1);

    data_nf2ff.Ephi = reshape(Ephi,s,1);
    data_nf2ff.E = [reshape(Ex,s,1),reshape(Ey,s,1),reshape(Ez,s,1)];
    data_nf2ff.Eabs = vecnorm(data_nf2ff.E,2,2);
    data_nf2ff.phi = reshape(phi_grid,s,1);
    data_nf2ff.theta= reshape(theta_grid,s,1);
    % ---------------------------------------------------------------------
end
