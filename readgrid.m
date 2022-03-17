function [x,nx]= readgrid(filename)
   
   fid = fopen(filename,'r');
    data  = fscanf(fid, '%f',[2 inf]);
    
    nx = data(1,1) ;
    x =data(1,2:end)'; 
   



end