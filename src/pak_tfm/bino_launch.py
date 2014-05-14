import b_sim
#b_sim.bino_simulation('inputdata_2_muchos.txt',30000)
displayinic=0
ciclos=10000
inputfile='4pim_ext'
haymut=0
haypred=0
haysup=0
el=0
plants_extinction={}
pols_extinction={} 
plants_extinction={'period':365,'spike':1,'start':0.05,'rate':-0.25,'numperiod':12,'species':[0,1,2]}
#pols_extinction={'period':365,'spike':0.7,'start':0.35,'rate':-0.25,'numperiod':3,'species':[0,1]}
pols_extinction={'period':365,'spike':1,'start':0.8,'rate':-0.15,'numperiod':100,'species':[2,0]}
dirsal='output/'
dirent='input/'
data_save=1
Na,Nb,Nc,Ra,Rb,maxa,maxb,maxreq,minreq=b_sim.bino_mutual(inputfile,ciclos,haymut,haypred,haysup,data_save,dirent,dirsal,eliminarenlaces=el,pl_ext=plants_extinction,pol_ext=pols_extinction)
if haypred:
    b_sim.food_render(Na,Nb,Nc,inputfile,displayinic,ciclos,dirsal,'')
b_sim.mutual_render(Na,Nb,Ra,Rb,maxa,maxb,maxreq,minreq,inputfile,displayinic,ciclos,dirsal,'')
#b_sim.phase_map(Na,Nb,1,1,displayinic,ciclos)plants_extinction