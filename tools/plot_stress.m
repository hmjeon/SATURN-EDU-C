clear all;

% Open file for post processing, fid is file pointer
fid = fopen('planestress_pos.txt','r');
n = fscanf(fid,'%i',1);
scalefactor = fscanf(fid,'%g',1);

% Setting the xy ratio of the figure, 1:1
set(gca,'DataAspectRatio',[1 1 1])

% Set up a color bar next to the figure
hc=colorbar;

% Use specific color map
colormap(jet(15));

% This holds the figure instead of creating new ones
% Useful to create overlapping 
hold on;

maxvalue = 0.0;
minvalue = 1.0e+14;

for j=1:n
   for i=1:4
      x(i) = fscanf(fid,'%g',1);
      y(i) = fscanf(fid,'%g',1);
      
      u(i) = fscanf(fid,'%g',1);
      v(i) = fscanf(fid,'%g',1);
      
      x(i) = x(i) + u(i) * scalefactor;
      y(i) = y(i) + v(i) * scalefactor;
      
      sxx(i) = fscanf(fid,'%g',1);
      syy(i) = fscanf(fid,'%g',1);
      sxy(i) = fscanf(fid,'%g',1);
      svm(i) = fscanf(fid,'%g',1);
   end

   xx = [x(1) x(2) x(3) x(4)]';
   yy = [y(1) y(2) y(3) y(4)]';
   vv = [svm(1) svm(2) svm(3) svm(4)]';
   patch(xx,yy,vv, 'LineWidth',1.0);
   %patch(xx,yy,zz,vv, 'Marker','o','EdgeColor','g','LineWidth',1);
   %patch(xx,yy,zz,vv,'Marker','o','EdgeColor','black'...
   %     ,'LineWidth',3,'FaceColor','none');
   
   if(maxvalue < max(vv))
       maxvalue = max(vv);
   end
   
   if(minvalue > min(vv))
       minvalue = min(vv);
   end
  
end

maxvalue
minvalue

hold off


   