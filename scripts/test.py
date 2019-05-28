import ffip

dx = 5
dt = 5 * 0.5
m = ffip.FePt(frequency=1/800)


print(m.get_dis_epsilon(frequency=1/800, dt=dt))