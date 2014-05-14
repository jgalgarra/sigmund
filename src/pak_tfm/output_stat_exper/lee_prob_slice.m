function msal = lee_prob_slice(nomfich,cmin,cmax,sentido)
    Minput=dlmread(nomfich);
    Minput(:,1)=int32(Minput(:,1));
    Minput(:,2)=int32(Minput(:,2));
    ycoord = Minput(1,2,1);
    xsup = max(Minput(:,1));
    xinf = min(Minput(:,1));

    if (sentido=='H')
        for j=xsup+1:cmax
            fila = [j,ycoord,1.0];
            Minput = vertcat(Minput,fila);
        end

        for j=cmin:xinf-1
            fila = [j,ycoord,0.0];
            Minput = vertcat(fila,Minput);
        end
    else
        xcoord = Minput(1,1,1);
        ysup = max(Minput(:,2));
        yinf = min(Minput(:,2));
        ymax= cmax;
        ymin = cmin;
        for j=ysup+1:ymax
            fila = [xcoord,j,1.0];
            Minput = vertcat(Minput,fila);
        end

        for j=ymin:yinf-1
            fila = [xcoord,j,0.0];
            Minput = vertcat(fila,Minput);
        end
    end
    msal = Minput;
end