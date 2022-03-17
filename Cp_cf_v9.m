clear; format long; clc; %close all;
%% Input------------------------------------------------------------------
% Uniform body with more points in the nose

 filename = 'Current_pp.stl'; % PATH TO stl 
 dir = 'files/'; % DIR WITH TNO FILES
% 
 a = 2;   % Body semiaxis
 nx = 50; % Number of points in where x is computed
 nstart = 11;
 nend   = 99;
 stride = 1;
 Re =1/0.00166;
 rho = 1;
 U_inf = 0.05;

 invalid_pts = 0;
%% Read stl ---------------------------------------------------------------

nit = 1+(nend-nstart)/stride;

[p,t,tnorm] = stlread(filename,1); % (vertices of each point) , (point num for each point), ...
                                      %normal vector

                                   
trin=length(tnorm); % number of triangles
%p = p(invalid_pts+1:end,:);
%t(:,:) = t(:,:) - invalid_pts ;

normx=tnorm(:,1);normy=tnorm(:,2);normz=tnorm(:,3);

vertexx = p(:,1);
vertexy = p(:,2);
vertexz = p(:,3);
vertexc = zeros(trin,3); tri_area= zeros(trin,1); 
ignore = false(trin,1) ;
num = 0;
for i = 1:trin
    for j=1:3
        if (t(i,j)-invalid_pts <= 0), ignore(i) = true; num=num+1;break; end
    end
    if ~ignore(i)
        vertexc(i,:) = 1/3*(p(t(i,1),:)+p(t(i,2),:)+p(t(i,3),:));
        tri_area(i,1) = area3D([p(t(i,1),1) p(t(i,2),1) p(t(i,3),1)], ...
            [p(t(i,1),2) p(t(i,2),2) p(t(i,3),2)], ...
            [p(t(i,1),3) p(t(i,2),3) p(t(i,3),3)] ) ;
    end
end
disp(['There are ', num2str(num),' invalid triangles']);
disp(['So only ', num2str(trin-num), ' number of triangles will be used :('])
vertexcx = vertexc(:,1); vertexcy = vertexc(:,2); vertexcz = vertexc(:,3);
%scatter3(vertexcx,vertexcy,vertexcz) ;

%% Read pressure ----------------------------------------------------------
 
pres_avg=zeros(trin,1); pavg_total =0;

for n=nstart:stride:nend

    file     = [dir, sprintf('pres.imb0001.00%d.tno', n)] ;
    [p_pres,pres] = tnoread(file);
    
    for i = 1: length(p_pres)
      if ~ignore(p_pres(i))  
      pres_avg(p_pres(i)) =  pres(i)*normx(p_pres(i));
      
      end
    end  
end
% time averaged pressure

pres_avg = pres_avg/nit;


%% Compute Cp -----------------------------------------------------------

Cp =  (pres_avg)/(0.5*rho*(U_inf^2)) ;
cptot = 0; area_tot =0 ;

data_write = zeros(trin-num,4);
k=1;
for i=1:trin
    if ignore(i) || vertexcx(i)<0, Cp(i)=0; 
    else
    data_write(k,1) = vertexcx(i);
    data_write(k,2) = vertexcy(i);
    data_write(k,3) = vertexcz(i);
    data_write(k,4) = Cp(i);
    
    cptot = cptot + abs(Cp(i)*tri_area(i)*normx(i)) ;
    area_tot = area_tot+ abs(tri_area(i)*normx(i)) ;
    
    k=k+1;
    end
end

save('Cp_3D_Fr02','data_write')
%dlmwrite('Cp_3D.dat',data_write, 'delimiter','\t')
clear data_write
set(gcf, 'Position', [100, 100, 700, 500])
scatter3(vertexcx,vertexcy,vertexcz,2,Cp,'filled') ;
xlabel('Z'), ylabel('Y'), zlabel('X')
h = gca; 
h.FontSize = 15; colorbar; colormap jet 
view(90,0)  % XY
%disp(['Cp for the body is ', num2str(sum(abs(Cp))/k)]);
disp(['Cp for the body is ', num2str(cptot/area_tot)]);
print('Cp_scatter','-dpng','-r0') 
%% Compute Cf
% dwdz ------------------------------------------------------------
 dwdz_avg(1:trin) = 0.0; 
for n=nstart:stride:nend
%     read dwdz value   
count_nan=0;
    file = [dir, sprintf('dwdz.imb0001.00%d.tno', n)];
    [p_visc,dwdz]=tnoread(file);
    
    for i = 1: length(p_visc)
      if ~ignore(p_visc(i)) &&  ~isnan(dwdz(i)) 
        dwdz_avg(p_visc(i)) = dwdz_avg(p_visc(i)) + dwdz(i);
        %if isnan(dwdz(i)), count_nan=count_nan+1 ; end
      else
          count_nan=count_nan+1 ;
      end
    end
end
disp(['There are ', num2str(count_nan),' Nan values at valid points of dwdz']);
dwdz_avg = dwdz_avg/nit;

%----------------------------------------------------------------------
% dwdx
 dwdx_avg(1:trin) = 0.0; 
for n=nstart:stride:nend
   count_nan=0;
%     read dwdx value    
    file = [dir, sprintf('dwdx.imb0001.00%d.tno', n)];
    [p_visc,dwdx]=tnoread(file);
    
    for i = 1: length(p_visc)
      if ~ignore(p_visc(i))  && ~isnan(dwdx(i))
        dwdx_avg(p_visc(i)) = dwdx_avg(p_visc(i)) + dwdx(i);
      else
          count_nan=count_nan+1 ;
      end
    end
end
disp(['There are ', num2str(count_nan),' Nan values at valid points of dwdx']);
dwdx_avg = dwdx_avg/nit;

%----------------------------------------------------------------------
% dwdy
 dwdy_avg(1:trin) = 0.0; count_nan=0;
for n=nstart:stride:nend
count_nan=0;
%     read dwdx value (always positive)    
    file = [dir, sprintf('dwdy.imb0001.00%d.tno', n)];
    [p_visc,dwdy]=tnoread(file);
    
    for i = 1: length(p_visc)
      if ~ignore(p_visc(i))  && ~isnan(dwdy(i))
        dwdy_avg(p_visc(i)) = dwdy_avg(p_visc(i)) + dwdy(i);
      else
          count_nan=count_nan+1 ;
      end
    end
end
disp(['There are ', num2str(count_nan),' Nan values at valid points of dwdy']);
dwdy_avg = dwdy_avg/nit;

%% ----------------------------------------------------------------------
% Donot run this section unless u need all components of drag
%{
% dudx
 dudx_avg(1:trin) = 0.0; count_nan=0;
for n=nstart:stride:nend
%     Store sign of dudx
    file = [dir, sprintf('sign_dudx.imb0001.000%d.tno', n)];
    [p_visc,sign_dudx]=tnoread(file);
    
    for i = 1:trin
        if sign_dudx(i,1) ~= 0
            sign_dudx(i,1)=sign_dudx(i,1)/abs(sign_dudx(i,1));
        else
            sign_dudx(i,1) = 0;
        end
    end
    
%     read dudx value (always positive)    
    file = [dir, sprintf('dudx.imb0001.000%d.tno', n)];
    [~,dudx]=tnoread(file);
    
    for i = 1: trin
      if ~ignore(p_visc(i))  
        dudx_avg(p_visc(i)) = dudx_avg(p_visc(i)) + dudx(i)*sign_dudx(i);
        if isnan(dudx(i)), count_nan=count_nan+1 ; end
      end
    end
end
disp(['There are ', num2str(count_nan),' Nan values at valid points of dudx']);
dudx_avg = dudx_avg/nit;

%----------------------------------------------------------------------
% dudy
 dudy_avg(1:trin) = 0.0; count_nan=0;
for n=nstart:stride:nend
%     Store sign of dudy
    file = [dir, sprintf('sign_dudy.imb0001.000%d.tno', n)];
    [p_visc,sign_dudy]=tnoread(file);
    
    for i = 1:trin
        if sign_dudy(i,1) ~= 0
            sign_dudy(i,1)=sign_dudy(i,1)/abs(sign_dudy(i,1));
        else
            sign_dudy(i,1) = 0;
        end
    end
    
%     read dudx value (always positive)    
    file = [dir, sprintf('dudy.imb0001.000%d.tno', n)];
    [~,dudy]=tnoread(file);
    
    for i = 1: trin
      if ~ignore(p_visc(i))  
        dudy_avg(p_visc(i)) = dudy_avg(p_visc(i)) + dudy(i)*sign_dudy(i);
        if isnan(dudy(i)), count_nan=count_nan+1 ; end
      end
    end
end
disp(['There are ', num2str(count_nan),' Nan values at valid points of dudy']);
dudy_avg = dudy_avg/nit;

%----------------------------------------------------------------------
% dvdy
 dvdy_avg(1:trin) = 0.0; count_nan=0;
for n=nstart:stride:nend
%     Store sign of dvdy
    file = [dir, sprintf('sign_dvdy.imb0001.000%d.tno', n)];
    [p_visc,sign_dvdy]=tnoread(file);
    
    for i = 1:trin
        if sign_dvdy(i,1) ~= 0
            sign_dvdy(i,1)=sign_dvdy(i,1)/abs(sign_dvdy(i,1));
        else
            sign_dvdy(i,1) = 0;
        end
    end
    
%     read dvdy value (always positive)    
    file = [dir, sprintf('dvdy.imb0001.000%d.tno', n)];
    [~,dvdy]=tnoread(file);
    
    for i = 1: trin
      if ~ignore(p_visc(i))  
        dvdy_avg(p_visc(i)) = dvdy_avg(p_visc(i)) + dvdy(i)*sign_dudy(i);
        if isnan(dvdy(i)), count_nan=count_nan+1 ; end
      end
    end
end
disp(['There are ', num2str(count_nan),' Nan values at valid points of dvdy']);
dvdy_avg = dvdy_avg/nit;

%}
% ------------------------------------------------------------------
%% Compute Cf -----------------------------------------------------------

Cf(1:trin)=0;
data_write = zeros(trin-num,4);
k=1; cftot = 0; 

for i=1:trin
    if ignore(i) || vertexcx(i)<0, Cf(i)=0; 
    else
    data_write(k,1) = vertexcx(i);
    data_write(k,2) = vertexcy(i);
    data_write(k,3) = vertexcz(i);
    data_write(k,4) = (dwdx_avg(i)*normz(i) + dwdy_avg(i)*normy(i) + ...
        dwdz_avg(i)*normx(i));
    data_write(k,4) = data_write(k,4)*4/(rho*U_inf^2)  ;
    
    Cf(i) =   data_write(k,4);
    cftot = cftot + abs(Cf(i)*tri_area(i)*normx(i)) ;
        
    k=k+1;
    end
end
set(gcf, 'Position', [100, 100, 700, 500])
scatter3(vertexcx,vertexcy,vertexcz,2,Cf/Re,'filled') ;
xlabel('Z'), ylabel('Y'), zlabel('X')
h = gca;
h.FontSize = 15;colorbar; colormap jet 
view(90,0)  % XY
%disp(['Cf for the body is ', num2str(sum(abs(Cf))/k)]);
disp(['Cf for the body is ', num2str(cftot/area_tot)]);
save('Cf_3D_Fr02','data_write')

print('Cf_scatter','-dpng','-r0')
%dlmwrite('Cf_3D_Fr05.dat',data_write, 'delimiter','\t')
