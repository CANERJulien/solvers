Tsatu=scalar(1)/(scalar(1)/T0-(scalar(8.31447)/(M*lv))*(log(p)-log(P0)));

mc=scalar(2)*fc/(scalar(2)-fc)*sqrt(M/(scalar(2*3.141459265358932384*8.31447)*T0))*P0*M*lv/(scalar(8.31447)*sqr(T0))*max(Tsatu-T,scalar(0));


mv=scalar(2)*fv/(scalar(2)-fv)*sqrt(M/(scalar(2*3.141459265358932384*8.31447)*T0))*P0*M*lv/(scalar(8.31447)*sqr(T0))*min(Tsatu-T,scalar(0));

m=mv+mc;

