{
//pSat = P0*exp((molWeight*deltaHvap/R)*(1/T0-1/mixture.T()));
//pSat = P0*exp((mixture.thermo2().W()*deltaHvap/R)*(1/T0-1/mixture.T()));
//pSat = P0*exp((W*deltaHvap/R)*(1/T0-1/mixture.T()));
pSat = P0*exp((W*deltaHvap/R)*(1/T0-1/T));
}
