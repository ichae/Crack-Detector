function showCurveAndPhi(phi,Img,cl)
	imshow(Img,[],'InitialMagnification',100);
    hold on; 
	[c,h] = contour(phi,[0 0],cl,'Linewidth',3); hold off;
	test = isequal(size(c,2),0);
	while (test==false)
        s = c(2,1);
        if ( s == (size(c,2)-1) )
            t = c;
            hold on; plot(t(1,2:end)',t(2,2:end)',cl,'Linewidth',3);
            test = true;
        else
            t = c(:,2:s+1);
            hold on; plot(t(1,1:end)',t(2,1:end)',cl,'Linewidth',3);
            c = c(:,s+2:end);
        end
	end    



end

