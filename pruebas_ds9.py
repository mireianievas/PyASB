import pyfits
import ds9
import numpy as np


def cargar_fits(fichero):
	print 'Cargando fichero fits'
	ficherofits = pyfits.open(fichero)
	datosfits = np.array(ficherofits[0].data, dtype='uint16')
	cabecerafits = ficherofits[0].header
	return datosfits,cabecerafits

def sesion_ds9(npdata):
	print 'Creando sesion ds9'
	myds9 = ds9.ds9()
	myds9.set_np2arr(npdata)
	myds9.set('zscale')
	myds9.set('zoom to fit')
	myds9.set('smooth radius 2')
	myds9.set('smooth yes')
	return myds9

def leer_coordenadas(myds9):
	print 'Pincha en una estrella en ds9'
	x,y = map(float, myds9.get('imexam any coordinate image').split(" ")[-2:])
	return x,y





#if __name__ == '__main__':
datosfits,cabecerafits = cargar_fits("/home/minaya/Johnson_V20120818_043106.fit")
myds9 = sesion_ds9(npdata=datosfits)

leer_coord = 'yes'
x = []; y = []
while leer_coord == 'yes':
	leer_coord = raw_input('Quiere leer una nueva coordenada? [yes,no]:')
	if leer_coord == 'yes':
		x_,y_ = leer_coordenadas(myds9)
		x.append(x_)
		y.append(y_)

print x,y

