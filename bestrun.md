gammafirst3nostack2
stacking='nostacking'
min_nucleation=1
MOD.setgeometry(theta=120, phi = 250)

OPT = Options(Nsim=3000)

H = HDNA(data, EXPNAME, model=MOD, options=OPT)
H.run([6e7, 3e7])

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>