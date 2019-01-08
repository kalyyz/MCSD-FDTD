function show(nx,ny,ez)
for i=1:nx 
 for j=1:ny 

   field(i,j)=real(ez(i,j,3)) ; 

 end ; 
end ;

mesh(1:nx,1:ny, real(field'))

axis([1 80 1 100 -5 5])
drawnow;
end