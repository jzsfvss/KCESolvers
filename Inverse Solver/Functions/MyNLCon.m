function [ ef, eqc ] = MyNLCon(EF, ku)

ef = EF([ ku(1), ku(2) ]) - ku(3);
eqc = 0;