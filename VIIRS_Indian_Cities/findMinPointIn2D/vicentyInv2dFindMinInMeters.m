% MatLab code for Vicenty's formula, which can be found at the link below: 
% https://en.wikipedia.org/wiki/Vincenty%27s_formulae
function [I, J, vicenty_ch, iterations] = vicentyInv2dFindMinInMeters(Lat, Lon, lat, lon, maxIter)
    a=6378137.0;
    f=1/298.257223563;
    b=(1-f).*a;
    
    phi_1 = Lat; L_1 = Lon;
    phi_2 = lat; L_2 = lon; 

    u_1= atan((1-f).*tan(deg2rad(phi_1)));
    u_2= atan((1-f).*tan(deg2rad(phi_2)));

    L=deg2rad(L_2-L_1);

    Lambda=L;  

    sin_u1=sin(u_1);
    cos_u1=cos(u_1);
    sin_u2=sin(u_2);
    cos_u2=cos(u_2);

    for cnt = 1 : maxIter
        cos_lambda=cos(Lambda);
        sin_lambda=sin(Lambda);
        sin_sigma=sqrt((cos_u2.*sin(Lambda)).^2+(cos_u1.*sin_u2-sin_u1.*cos_u2.*cos_lambda).^2);
        cos_sigma=sin_u1.*sin_u2+cos_u1.*cos_u2.*cos_lambda;
        sigma=atan2(sin_sigma,cos_sigma);
        sin_alpha=(cos_u1.*cos_u2.*sin_lambda)./sin_sigma;
        cos_sq_alpha=1-sin_alpha.^2;
        cos2_sigma_m=cos_sigma-((2.*sin_u1.*sin_u2)./cos_sq_alpha);
        C=(f./16).*cos_sq_alpha.*(4+f.*(4-3.*cos_sq_alpha));
        Lambda_prev=Lambda;
        Lambda=L+(1-C).*f.*sin_alpha.*(sigma+C.*sin_sigma.*(cos2_sigma_m+C.*cos_sigma.*(-1+2.*cos2_sigma_m.^2)));

        diff=abs(Lambda_prev-Lambda);
        if all(diff(:))<=10^-12
            break;
        end

        if cnt == maxIter
            error(['Failed to converge in ', num2str(maxIter)]);
        end
    end

    u_sq=cos_sq_alpha.*((a^2-b^2)./b^2);
    A=1+(u_sq./16384).*(4096+u_sq.*(-768+u_sq.*(320-175.*u_sq)));
    B=(u_sq./1024).*(256+u_sq.*(-128+u_sq.*(74-47.*u_sq)));
    delta_sig=B.*sin_sigma.*(cos2_sigma_m+0.25.*B.*(cos_sigma.*(-1+2.*cos2_sigma_m.^2)-(1./6).*B.*cos2_sigma_m.*(-3+4.*sin_sigma.^2).*(-3+4.*cos2_sigma_m.^2)));
    distance = b.*A.*(sigma-delta_sig);
    [I,J]=find(distance(:,:)==min(distance(:))); 
    vicenty_ch = distance(I,J);
    iterations = cnt;
end