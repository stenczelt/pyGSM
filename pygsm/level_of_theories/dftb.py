# standard library imports
import os
import re
import subprocess

# third party
import numpy as np

# local application imports
from pygsm import utilities
from .base_lot import Lot


class DFTB(Lot):

    def __init__(self,options):
        super(DFTB,self).__init__(options)
        os.system('rm -f dftb_jobs.txt')
        print(" making folder scratch/{}".format(self.node_id))
        os.system('mkdir -p scratch/{}'.format(self.node_id))
        os.system('cp {} scratch/{}'.format(self.lot_inp_file,self.node_id))

    def run(self,geom):
        owd = os.getcwd()
        utilities.utilities.manage_xyz.write_xyz('scratch/{}/tmp.xyz'.format(self.node_id), geom, scale=1.0)
        os.system('./xyz2gen scratch/{}/tmp.xyz'.format(self.node_id))
        os.chdir('scratch/{}'.format(self.node_id))
        os.system('pwd')
        cmd = "dftb+"
        proc = subprocess.Popen(cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                )
        stdout,stderr = proc.communicate()
        #with open('dftb_jobs.txt','a') as out:
        #    out.write(stdout)
        #    out.write(stderr)

        ofilepath = "detailed.out"
        with open(ofilepath,'r') as ofile:
            olines = ofile.readlines()

        self.E = []
        temp = 0
        tmpgrada=[]
        tmpgrad=[]
        pattern=re.compile(r"Total energy:                     [-+]?[0-9]*\.?[0-9]+ H")
        for line in olines:
            for match in re.finditer(pattern,line):
                tmpline = line.split()
                self.E.append((1,0,float(tmpline[2])))
            if line==" Total Forces\n":
                temp+=1
            elif temp>0:
                tmpline = line.split()
                tmpgrad.append([float(i) for i in tmpline])
                temp+=1
            if temp> len(self.atoms):
                break
        tmpgrada.append(tmpgrad)

        self.grada=[]
        for count,i in enumerate(self.states):
            if i[0]==1:
                self.grada.append((1,i[1],tmpgrada[count]))
            if i[0]==3:
                self.grada.append((3,i[1],tmpgrada[count]))
        self.hasRanForCurrentCoords=True

        os.chdir(owd)
        return

    def get_energy(self,coords,multiplicity,state):
        if self.hasRanForCurrentCoords==False or (coords != self.currentCoords).any():
            self.currentCoords = coords.copy()
            geom = utilities.utilities.manage_xyz.np_to_xyz(self.geom, self.currentCoords)
            self.run(geom)
        tmp = self.search_PES_tuple(self.E,multiplicity,state)[0][2]
        return self.search_PES_tuple(self.E,multiplicity,state)[0][2] * utilities.utilities.units.KCAL_MOL_PER_AU

    def get_gradient(self,coords,multiplicity,state):
        if self.hasRanForCurrentCoords==False or (coords != self.currentCoords).any():
            self.currentCoords = coords.copy()
            geom = utilities.utilities.manage_xyz.np_to_xyz(self.geom, self.currentCoords)
            self.run(geom)
        tmp = self.search_PES_tuple(self.grada,multiplicity,state)[0][2]
        return np.asarray(tmp) * -1. * utilities.utilities.units.ANGSTROM_TO_AU

if __name__=='__main__':
    filepath="../../data/ethylene.xyz"
    dftb = DFTB.from_options(states=[(1,0)],fnm=filepath,lot_inp_file='../../data/dftb_in.hsd')
    geom=utilities.manage_xyz.read_xyz(filepath)
    xyz = utilities.manage_xyz.xyz_to_np(geom)
    print(dftb.get_energy(xyz,1,0))
    print(dftb.get_gradient(xyz,1,0))

    xyz = xyz+ np.random.rand(xyz.shape[0],xyz.shape[1])*0.1
    print(dftb.get_energy(xyz,1,0))
    print(dftb.get_gradient(xyz,1,0))
