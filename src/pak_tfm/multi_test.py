from multiprocessing import Pool

class DummyC():
    def __init__(self,a,b):
        self.a = a
        self.b = b


def f(x):
    return x.b*2+x.a*3

if __name__ == '__main__':
    val=[]
    for i in range(0,3):
        ins = DummyC(i,2*i)
        val.append(ins)
    with Pool(5) as p:
        print(p.map(f, val))
