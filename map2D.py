import numpy as np

class map2D(object):
    def __init__(self, twiss=[1.0,0.0], tune=0.0, chrom=0.0,
                 npart=1, emit=1e-6, twiss_beam=None, espr=0.0, particles=None):

        self.set_map(twiss, tune, chrom)
        self.set_particle(npart, emit, twiss_beam, espr, particles)

    def set_map(self, twiss, tune, chrom):
        '''
        Set the 2-by-2 map.
        :param twiss: [beta, alpha] in a list
        :param tune: The tune of map
        :param chrom: set chromaticity
        :return:
        '''
        self.beta=twiss[0]
        self.alpha=twiss[1]
        self.gamma=(1+self.alpha*self.alpha)/self.beta
        self.chrom=chrom
        self.phi=2*np.pi*tune
        self.linmap=np.array([[np.cos(self.phi)+self.alpha*np.sin(self.phi), self.beta*np.sin(self.phi)],
                              [-np.sin(self.phi)*self.gamma, np.cos(self.phi) - self.alpha * np.sin(self.phi)]])


    def set_particle(self, npart=1, emit=1e-6, twiss_beam=None, espr=0.0, particles=None):
        '''
        Set the particles either by external input or random generation
        :param np:  Number of particles default is 1, will be ignored if particles is given
        :param emit: rms emittance of the beam, default is 1e-6 m-rad
        :param twiss_beam:[beta, alpha] in a list.  if None, matched beam will be generated
        :param espr: energy spread, default is zero
        :param particles:  External given particle array in shate of (2,N)
        :return:
        '''
        if particles is None:
            self.npart=npart
            if twiss_beam is None:
                self.beam_beta=self.beta
                self.beam_alpha=self.alpha
            else:
                self.beam_beta = twiss_beam[0]
                self.beam_alpha = twiss_beam[1]
            self.coor2D=np.zeros((2, npart))
            self.coor2D[0, :]=np.random.randn(npart)
            avex=np.average(self.coor2D[0, :])
            sizex=np.std(self.coor2D[0, :])
            self.coor2D[0, :] -= avex
            self.coor2D[0, :] *= (np.sqrt(emit * self.beam_beta) / sizex)

            self.coor2D[1, :] = np.random.randn(npart)
            avexp = np.average(self.coor2D[1, :])
            sizexp = np.std(self.coor2D[1, :])
            self.coor2D[1, :] -= avexp
            self.coor2D[1, :] *= (np.sqrt(emit / self.beam_beta) / sizexp)

            cxxp=np.average((self.coor2D[0, :]) * (self.coor2D[1, :]))
            emit=np.sqrt(sizex * sizexp * sizex * sizexp - cxxp * cxxp)


            self.coor2D[1, :]+= self.coor2D[0, :] * ((-self.beam_alpha + cxxp / emit) / self.beam_beta)


        else:
            self.npart = particles.shape[1]
            self.coor2D=particles

        if espr!=0 and self.chrom!=0:
            self.espread=np.random.randn(self.npart)*espr
            phi = self.phi + (2 * np.pi * self.chrom) * self.espread
            csphi = np.cos(phi)
            snphi = np.sin(phi)
            self.chrom_m11 = csphi + self.alpha * snphi
            self.chrom_m12 = self.beta * snphi
            self.chrom_m21 = -self.gamma * snphi
            self.chrom_m22 = csphi - self.alpha * snphi
        else:
            self.espread=None


    def propagate(self):
        if self.espread is None:
            self.coor2D= self.linmap @ self.coor2D
        else:
            self.coor2D[0, :], self.coor2D[1, :] = self.chrom_m11 * self.coor2D[0, :] + self.chrom_m12 * self.coor2D[1, :], self.chrom_m21 * self.coor2D[0, :] + self.chrom_m22 * self.coor2D[1, :]


    def statistics(self):
         avex = np.average(self.coor2D[0, :])
         avexp = np.average(self.coor2D[1, :])
         sizex = np.std(self.coor2D[0, :])
         sizexp = np.std(self.coor2D[1, :])
         cxxp= np.average((self.coor2D[0, :] - avex) * (self.coor2D[1, :] - avexp))
         return avex, avexp, sizex, sizexp, np.sqrt(sizex*sizexp*sizex*sizexp -cxxp*cxxp)