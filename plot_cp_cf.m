% Make plots of Cp and Cf on a structured grid and plot them
% I am switching coordinate system here

clc ; clear ;

load('Cp_3D_Fr02.mat');
sorted = sortrows(data_write,1) ;
x = sorted(:,3) ;
y = sorted(:,2) ;
z = sorted(:,1) ;
cp = sorted(:,4) ;

dz =0.5;
nbeta = 256 ;
zmax = 150 ;
bool(1:3) = true;

null = 0; num = 0; zp = -2;
for k = 1:length(z)
  if  z(k) > zp + dz
      r = sqrt(x(k)^2 + y(k)^2) ;
        
      if z(k)> zmax/4 && bool(1),nbeta = nbeta/2 ;dz = 1; bool(1)=false;end
      if z(k)> zmax/2 && bool(2),nbeta = nbeta/2 ;dz = 1; bool(2)=false;end  
      if z(k)> 3*zmax/4 && bool(3),nbeta = nbeta ;dz =0.5; bool(3)=false;end
        
      for i=num+1:num+nbeta
        znew(i) = z(k);
        xnew(i) = r*cos(2*pi*(i-num)/nbeta) ;
        ynew(i) = r*sin(2*pi*(i-num)/nbeta) ;
      end
      num = num + nbeta ;  
      zp = z(k) ;
   end
end
znew = znew' ; ynew = ynew';xnew =xnew' ;  

cp_new = griddata(x,y,z,cp,xnew,ynew,znew) ;
scatter3(xnew,ynew,znew,3,cp_new,'filled');
xlabel('x'), ylabel('y'), zlabel('z')