import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from Veryfunlib import veryfunlib as vfl


#setto il seed
vfl.setseed(12)
betac=0.44

#mi dice quanta roba c'è
for j in range(7):
	L=j+1
	a=str(L)				#converto in stringa l'indice, ho scritto i file in questo modo
	rep, camp, beta0, dbeta = np.loadtxt('shape'+a+'0.dat', unpack=True, max_rows=1)
	rep=int(rep)
	camp=int(camp)

	#print(camp)

	#carico i dati
	ene_data=np.zeros((camp, rep))
	for i in range(rep):
		L_dat = camp*i
		ene_data[:,i] = np.loadtxt('ene_'+a+'0.dat', unpack=True, skiprows=L_dat, max_rows=camp)



	#avvio il generatore random, usare MT19937 al posto di PCG64 per il mersenne twister


	#faccio le medie e le varianze
	Emean=np.zeros(rep)
	Evar=np.zeros(rep)
	Emean=np.mean(ene_data, axis=0)
	Evar=L*L*np.var(ene_data, axis=0)


	#bootstrap
	E10=np.zeros(rep)
	onetimeplot=False


	#controllo che sia abbastanza piccolo
	bootlen=5
	expect=camp/(2**bootlen)
	#print('campionamenti: %d' %camp )
	#print('bin: %d' %bootlen)

	if (expect<camp/100):
		print('bin troppo grande')
		print(expect)
		
	else:
		for i in range(rep):
			E10[i]=vfl.bootstrap(ene_data[:,i], bootlen, est='variance', rep=80, R_print=onetimeplot)
			 #onetimeplot=False
			 
	#plot chi	 
	X=np.linspace(beta0, beta0 + rep*dbeta , rep)
	
	#esponente critico capacità termica: alpha = 0
	esp=0
	#Mvar=Mvar/(L**esp)
	#Merr=Merr/(L**esp)
	X=X-betac
	X*=L
	
	plt.errorbar(X,Evar,E10, linestyle='', markersize=5, marker='.',label='L=' + a+'0')





plt.title('Capacità termica attorno al punto critico')
plt.ylabel('Capacità termica')
plt.xlabel(' Beta ')
plt.grid(color = 'silver')
plt.legend()
plt.show()


