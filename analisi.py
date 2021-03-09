import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from Veryfunlib import veryfunlib as vfl


#setto il seed
vfl.setseed(12)


#temperatura critica
betac=0.44

#mi dice quanta roba c'è
for j in range(7):
	L=j+1
	a=str(L)				#converto in stringa l'indice, ho scritto i file in questo modo
	rep, camp, beta0, dbeta = np.loadtxt('shape'+ a+ '0.dat', unpack=True, max_rows=1)			#carico il file con la forma dei dati
	rep=int(rep)
	camp=int(camp)

	#print(camp)

	#carico i dati
	mag_data=np.zeros((camp, rep))
	for i in range(rep):														
		L_dat = camp*i
		mag_data[:,i] = np.loadtxt('mag_'+a+'0.dat', unpack=True, skiprows=L_dat, max_rows=camp)



	#avvio il generatore random, usare MT19937 al posto di PCG64 per il mersenne twister


	#faccio le medie e le varianze
	Mmean=np.zeros(rep)
	Mvar=np.zeros(rep)
	Mmean=np.mean(abs(mag_data), axis=0)
	Mvar=L*L*np.var(abs(mag_data), axis=0)

	#print('medie')
	#print(len(Mmean))
	#print(Mmean)
	#print('-')

	#bootstrap
	Merr=np.zeros(rep)
	onetimeplot=False
	for i in range(rep):
		Merr[i]=vfl.bootstrap(mag_data[:,i], 7, est='mean', rep=80, R_print=onetimeplot)
		onetimeplot=False
		
	X=np.linspace(beta0, beta0 + rep*dbeta , rep)
	#Merr=L*L*Merr
	#print(len(X))
	#print(len(M10))
	#print(len(Mvar))
	#print(M10)


	
	
	#finite size scaling gamma:7/4, beta:1/8 (entra con un - nella magnetizzazione)
	esp=-1/8
	Mvar=Mvar/(L**esp)
	Merr=Merr/(L**esp)
	X=X-betac
	X*=L
	
	#plot chi
	plt.errorbar(X, Mmean, Merr, linestyle='', markersize=5, marker='.', label='L=' + a+'0')
	
plt.title('Magnetizzazione attorno al punto critico')
plt.ylabel('M')
plt.xlabel(' Beta ')
plt.grid(color = 'silver')
plt.legend()
plt.show()







