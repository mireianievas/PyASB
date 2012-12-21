import numpy as np
import ds9
import pyfits

class StarSurface:
	def __init__(self,pyfits_image):
		self.myds9 = ds9.ds9()
		self.pyfits_data = pyfits.open(pyfits_image)[0].data
		self.pyfits_filtered = pyfits.open(pyfits_image)[0].data
		self.starsurface = np.zeros(self.pyfits_data.shape,dtype='bool')
		self.myds9.set_np2arr(np.array(self.pyfits_data,dtype='uint16'))
		''' Seed for initiate iteration'''
		self.estimate_background()
		self.estimate_centroid()
		assert self.pyfits_data[self.ycentroid][self.xcentroid] > self.background
		self.starsurface[self.ycentroid,self.xcentroid] = True
		print 'Iterating: ',
		self.loop()
	
	def estimate_centroid(self):
		exponent = 2 # Valures > 1 intensify the convergence of the method.
		detection_threshold = 0.5
		pyfits_absdata = self.pyfits_data-self.background
		try:
			xweight = np.array([range(len(pyfits_absdata[0]))]*len(pyfits_absdata))
			yweight = np.array([range(len(pyfits_absdata))]*len(pyfits_absdata[0])).transpose()
			self.xcentroid = np.sum(xweight*pyfits_absdata**exponent)/np.sum(pyfits_absdata**exponent)
			self.ycentroid = np.sum(yweight*pyfits_absdata**exponent)/np.sum(pyfits_absdata**exponent)
			assert np.std(pyfits_absdata)>detection_threshold*np.mean(np.abs(pyfits_absdata))
		except:
			raise
	
	def estimate_background(self):
		self.background = np.median(self.pyfits_data[self.starsurface==False])
		#print 'Background: '+str(self.background)
	
	def build_star_test_pixels(self):
		self.testsurface = np.zeros(self.pyfits_data.shape,dtype='bool')
		for y in xrange(1,self.starsurface.shape[0]-1):
			for x in xrange(1,self.starsurface.shape[1]-1):
				if self.starsurface[y][x] == True:
					#print y,x
					if self.starsurface[y-1][x-1] == False: 
						self.testsurface[y-1,x-1] = True
					if self.starsurface[y+0][x-1] == False: 
						self.testsurface[y+0,x-1] = True
					if self.starsurface[y+1][x-1] == False: 
						self.testsurface[y+1,x-1] = True
					if self.starsurface[y-1][x+0] == False: 
						self.testsurface[y-1,x+0] = True
					if self.starsurface[y+1][x+0] == False: 
						self.testsurface[y+1,x+0] = True
					if self.starsurface[y-1][x+1] == False: 
						self.testsurface[y-1,x+1] = True
					if self.starsurface[y+0][x+1] == False: 
						self.testsurface[y+0,x+1] = True
					if self.starsurface[y+1][x+1] == False: 
						self.testsurface[y+1,x+1] = True
	
	def flux_increase_test(self):
		return np.mean(self.pyfits_data[self.testsurface==True]) > 1.01*self.background
	
	def incorporate_pixels(self):
		for y in xrange(self.starsurface.shape[0]):
			for x in xrange(self.starsurface.shape[1]):
				if self.testsurface[y][x]==True and self.pyfits_data[y][x]>self.background:
					self.starsurface[y][x]=True
				#if np.sum(self.starsurface[y-1:y+1,x-1:x+1])>2:
				#	self.starsurface[y][x]=True				
	
	def loop(self):
		print '.',
		self.pyfits_filtered[self.starsurface==True]=0
		#self.myds9.set_np2arr(np.array(self.pyfits_filtered,dtype='uint16'))
		self.estimate_background()
		self.build_star_test_pixels()
		#print self.testsurface
		if self.flux_increase_test() == True:
			self.incorporate_pixels()
			#print self.starsurface
			self.loop()
		
Star = StarSurface("recorte_estrella2.fits")
Star.myds9.set_np2arr(np.array(Star.pyfits_filtered,dtype='uint16'))
		
						
	
	