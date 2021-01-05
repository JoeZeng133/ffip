import ffip

ft = 1/500

m1 = ffip.Au
er1 = m1.get_dis_epsilon(ft, 1)

m2 = ffip.Au_IIR1
er2 = m2.get_dis_epsilon(ft, 1)

print(er1, er2)