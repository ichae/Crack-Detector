
function angle=givePrincipalArgumentAngle(x,y)

  if x==0 && y>0,
     angle=pi/2;
  elseif x==0 && y<0,
     angle=(3*pi)/2;
  elseif x>0 && y==0,
     angle=0;
  elseif x<0 && y==0,
     angle=pi;
  elseif x>0 && y>0, % first quadrant
      angle=atan(y/x);
  elseif x>0 && y<0, % fourth quadrant
      angle=2*pi-atan(abs(y)/x);
  elseif x<0 && y>0, % second quadrant
      angle=pi-atan(y/abs(x));
  else %  third quadrant
      angle=pi+atan(abs(y)/abs(x));
  end   
  
  return