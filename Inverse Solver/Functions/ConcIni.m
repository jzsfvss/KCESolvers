function ueq = ConcIni(Lini, Tini, Kd)

ueq = [ 0, 0, 0 ]';

a = Lini - Tini + (Kd/1000);
ueq(2) = 1000*0.5*abs(-a + sqrt((a^2) + Kd*Tini/250));
%ueq(2) = 1000*((-(Lini-Tini+Kd))+(sqrt(((Lini-Tini+Kd)^2)-(4*(-Kd*Tini)))))/2;
ueq(3) = 1000*Tini - ueq(2);
ueq(1) = 1000*Lini - ueq(3);