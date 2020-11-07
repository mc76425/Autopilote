import pandas as pa
import scipy.interpolate as sip

L = ["a0","Cx0","Cza","f"]

for i in L:
    df = pa.read_csv(i+".csv",names=["x","y"],delimiter=";")
    x, y = [float(j) for j in df.x.to_list()], [float(h) for h in df.y.to_list()]
    f = sip.interp1d(x, y)
    print(i," : ",f(1.4))