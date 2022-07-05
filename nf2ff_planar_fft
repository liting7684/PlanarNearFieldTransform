%%使用FFT实现平面波谱函数
  

function [data_nf2ff]=nf2ff_planar(Amp,M,N,lambda,dx,dy,deltax,deltay,Im,theta,phi,d,Mc,Nc,padding)
%构建结果表
p=length(theta);
t=length(phi);
data_nf2ff=table(zeros(p*t,1),zeros(p*t,1),zeros(p*t,3),zeros(p*t,1),zeros(p*t,1),zeros(p*t,1));
data_nf2ff.Properties.VariableNames={'theta','phi','E','Etheta','Ephi','Eabs'};
k0=2*pi/lambda;


%补零
Mc_padded=padding*Mc;
Nc_padded=padding*Nc;
m=-Mc_padded/2:Mc_padded/2-1;
n=-Nc_padded/2:Nc_padded/2-1;


 % --------------------------------------------------
 kx=2*pi*m/(Mc_padded*deltax);ky=2*pi*n/(Nc_padded*deltay);
 [ky_grid,kx_grid]=meshgrid(ky,kx);
 kz_grid=sqrt(k0^2-kx_grid.^2-ky_grid.^2);
 
  %球面远场波数矢量  套用公式直角坐标转换成球面坐标系
  [theta_grid,phi_grid]=meshgrid(theta,phi);%转换成网格向量  meshgrid 用来就是产生二维坐标矩阵
  kx_grid_spherical=k0*sin(theta_grid).*cos(phi_grid);
  ky_grid_spherical = k0*sin(theta_grid).*sin(phi_grid);
  kz_grid_spherical=k0*cos(theta_grid);
  
  
  % 先计算采样平面的Ex
  Ex=0;%初始化
  for m=0:Mc-1
      for n=0:Nc-1
          XL1=(m-Mc/2)*deltax;
          YL2=(n-Nc/2)*deltay;
          Ex(m+1,n+1)=Ex_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im);%计算Ex的采样电场值
      end
  end
%   surf(abs(Ex));%输出Ex的幅值和相位
%   surf(angle(Ex));
  
  %计算采样平面的Ey
  Ey=0;
    for m=0:Mc-1
        for n=0:Nc-1
            XL1=(m-Mc/2)*deltax;
            YL2=(n-Nc/2)*deltay;
            Ey(m+1,n+1)=Ey_nfield(Amp,M,N,lambda,dx,dy,XL1,YL2,d,Im,Mc,Nc,deltax,deltay);%计算Ey的采样电场值
        end
    end
%   surf(abs(Ey));%输出Ey的幅值和相位
%   surf(angle(Ey));
  %--------------------------------------------------------------
  Ax=ifftshift(ifft2(Ex,Mc_padded,Nc_padded));
  Ay=ifftshift(ifft2(Ey,Mc_padded,Nc_padded));
  Az=-(Ax.*kx_grid+Ay.*ky_grid)./kz_grid;
  surf(abs(Ax));
  figure(20)
  surf(abs(Ay));
  title('Ay');
  %--------------------------------------------------------------
  [kx,ky]=meshgrid(ky,kx);
  % Interpolate Modes in spherical coordinates
  Ax_ff=interp2(kx,ky,abs(Ax),kx_grid_spherical,ky_grid_spherical,'spline');
  Ay_ff=interp2(kx,ky,abs(Ay),kx_grid_spherical,ky_grid_spherical,'spline');
  Az_ff=interp2(kx,ky,abs(Az),kx_grid_spherical,ky_grid_spherical,'spline');
  
  
  
  % Far Field 
    r=100;
    C=1j*(k0*exp(-1j*k0*r))/(2*pi*r);
    Etheta=C*(Ax_ff.*cos(phi_grid)+Ay_ff.*sin(phi_grid));
    Ephi=C*cos(theta_grid).*(-Ax_ff.*sin(phi_grid)+Ay_ff.*cos(phi_grid));
    Ex = C*cos(theta_grid).*Ax_ff;
    Ey = C*cos(theta_grid).*Ay_ff*0;
    Ez = C*cos(theta_grid).*Az_ff*0;

    % --------------------------------------------------
    s = numel(theta_grid);
    data_nf2ff.Etheta = reshape(Etheta,s,1);
    data_nf2ff.Ephi = reshape(Ephi,s,1);
    data_nf2ff.E = [reshape(Ex,s,1),reshape(Ey,s,1),reshape(Ez,s,1)];
    data_nf2ff.Eabs = vecnorm(data_nf2ff.E,2,2);
    data_nf2ff.phi = reshape(phi_grid,s,1);
    data_nf2ff.theta= reshape(theta_grid,s,1);
    % --------------------------------------------------
end

  
  
  
  
  
  
  
  
  
  
  
 
