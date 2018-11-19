import numpy as np
import options
import os
from base_gsm import *
from dlc import *
from pes import *
import pybel as pb
import sys

class GSM(Base_Method):
    @staticmethod
    def from_options(**kwargs):
        return GSM(GSM.default_options().set_values(kwargs))

    def __init__(
            self,
            options,
            ):

        super(GSM,self).__init__(options)

        tmp = self.options['ICoord2']
        self.icoords[0] = DLC.union_ic(self.icoords[0],tmp)
        print "after union"
        lot1 = tmp.PES.lot.copy(
                tmp.PES.lot, 
                self.nnodes-1)
        PES1 = PES(tmp.PES.options.copy().set_values({
            "lot": lot1,
            }))
        self.icoords[-1] = DLC(self.icoords[0].options.copy().set_values(dict(
            mol= tmp.mol,
            PES=PES1,
            )))
        print "print levels at beginning are ",self.icoords[0].print_level
        print "print levels at beginning are ",self.icoords[-1].print_level

    def go_gsm(self,g_iters=1,maxopt=1,nconstraints=1,o_iters=1,optsteps=1):
        self.icoords[0].gradrms = 0.
        self.icoords[-1].gradrms = 0.
        self.icoords[0].energies = self.icoords[0].PES.get_energy(self.icoords[0].geom)
        self.icoords[-1].energy = self.icoords[-1].PES.get_energy(self.icoords[-1].geom)
        self.interpolate(2) 
        self.growth_iters(iters=g_iters,maxopt=maxopt,nconstraints=nconstraints)
        self.opt_iters(max_iter=o_iters,optsteps=optsteps)


    def interpolate(self,newnodes=1):
        if self.nn+newnodes > self.nnodes:
            print("Adding too many nodes, cannot interpolate")
        sign = -1
        for i in range(newnodes):
            print "Adding node",i
            sign *= -1
            if sign == 1:
                self.interpolateR()
            else:
                self.interpolateP()

    def add_node(self,n1,n2,n3):
        print "adding node: %i between %i %i" %(n2,n1,n3)
        return DLC.add_node(self.icoords[n1],self.icoords[n3],self.nnodes,self.nn)

    def set_active(self,nR,nP):
        print(" Here is active:",self.active)
        if nR!=nP:
            print(" setting active nodes to %i and %i"%(nR,nP))
        else:
            print(" setting active node to %i "%nR)

        for i in range(self.nnodes):
            if self.icoords[i] !=0:
                self.active[i] = False;
                self.icoords[i].OPTTHRESH = self.CONV_TOL*2.;
        self.active[nR] = True
        self.active[nP] = True
        print(" Here is new active:",self.active)

    def tangent(self,n1,n2):
        print" getting tangent from between %i %i pointing towards %i"%(n2,n1,n2)
        return DLC.tangent_1(self.icoords[n2],self.icoords[n1])

    def check_if_grown(self):
        isDone=True
        for act in self.active:
            if act:
                isDone = False
                break
        return isDone

    def check_add_node(self):
        if self.icoords[self.nR-1].gradrms < self.gaddmax:
            self.active[self.nR-1] = False
            if self.icoords[self.nR] == 0:
                self.interpolateR()
                print "print_level is",self.icoords[self.nR-1].print_level
        if self.icoords[self.nnodes-self.nP].gradrms < self.gaddmax:
            self.active[self.nnodes-self.nP] = False
            if self.icoords[-self.nP-1] == 0:
                self.interpolateP()
                print "print_level is",self.icoords[self.nR-1].print_level

    def make_nlist(self):
        ncurrent = 0
        nlist = [0]*(2*self.nnodes)
        for n in range(self.nR-1):
            nlist[2*ncurrent] = n
            nlist[2*ncurrent+1] = n+1
            ncurrent += 1

        for n in range(self.nnodes-self.nP+1,self.nnodes):
            nlist[2*ncurrent] = n
            nlist[2*ncurrent+1] = n-1
            ncurrent += 1

        nlist[2*ncurrent] = self.nR -1
        nlist[2*ncurrent+1] = self.nnodes - self.nP

        if False:
            nlist[2*ncurrent+1] = self.nR - 2 #for isMAP_SE

        #TODO is this actually used?
        if self.nR == 0: nlist[2*ncurrent] += 1
        if self.nP == 0: nlist[2*ncurrent+1] -= 1
        ncurrent += 1
        nlist[2*ncurrent] = self.nnodes -self.nP
        nlist[2*ncurrent+1] = self.nR-1
        #TODO is this actually used?
        if self.nR == 0: nlist[2*ncurrent+1] += 1
        if self.nP == 0: nlist[2*ncurrent] -= 1
        ncurrent += 1

        return ncurrent,nlist

    def start_string(self):
        print "\n"
        self.interpolate(2) 
        self.nn=2
        self.nR=1
        self.nP=1

if __name__ == '__main__':
#    from icoord import *
    ORCA=False
    QCHEM=True
    PYTC=False
    nproc=8

    if QCHEM:
        from qchem import *
    if ORCA:
        from orca import *
    if PYTC:
        from pytc import *
    import manage_xyz

    if False:
        filepath="tests/fluoroethene.xyz"
        filepath2="tests/stretched_fluoroethene.xyz"
        nocc=11
        nactive=2

    if True:
        filepath2="tests/SiH2H2.xyz"
        filepath="tests/SiH4.xyz"
        nocc=8
        nactive=2

    mol=pb.readfile("xyz",filepath).next()
    mol2=pb.readfile("xyz",filepath2).next()
    basis = 'sto-3g'
    if QCHEM:
        lot=QChem.from_options(states=[(1,0)],charge=0,basis=basis,functional='B3LYP',nproc=nproc)
        lot2=QChem.from_options(states=[(1,0)],charge=0,basis=basis,functional='B3LYP',nproc=nproc)
    
    if ORCA:
        lot=Orca.from_options(states=[(1,0)],charge=0,basis=basis,functional='wB97X-D3',nproc=nproc)
        lot2=Orca.from_options(states=[(1,0)],charge=0,basis=basis,functional='wB97X-D3',nproc=nproc)
    if PYTC:
        lot=PyTC.from_options(states=[(1,0)],nocc=nocc,nactive=nactive,basis='6-31gs')
        lot.cas_from_file(filepath)
        lot2=PyTC.from_options(states=[(1,0)],nocc=nocc,nactive=nactive,basis='6-31gs')
        lot2.casci_from_file_from_template(filepath,filepath2,nocc,nocc)

    pes = PES.from_options(lot=lot,ad_idx=0,multiplicity=1)
    pes2 = PES.from_options(lot=lot2,ad_idx=0,multiplicity=1)

    print "\n IC1 \n"
    ic1=DLC.from_options(mol=mol,PES=pes,print_level=0)
    print "\n IC2 \n"
    ic2=DLC.from_options(mol=mol2,PES=pes2,print_level=0)
        
    nnodes=9
    if True:
        print "\n Starting GSM \n"
        gsm=GSM.from_options(ICoord1=ic1,ICoord2=ic2,nnodes=nnodes,nconstraints=1,CONV_TOL=0.001)
        gsm.go_gsm(g_iters=30,maxopt=1,nconstraints=1,o_iters=30,optsteps=1)

    if False:
        print DLC.tangent_1(gsm.icoords[0],gsm.icoords[-1])
    
    if False:
        for i in range(gsm.nnodes):
            gsm.icoords[i].energy = gsm.icoords[i].PES.get_energy(gsm.icoords[i].geom)
        gsm.get_tangents_1e(n0=0)

    if False:
        gsm.ic_reparam_g()


    if False:
        gsm.grow_string(50)
        #gsm.growth_iters(iters=50,maxopt=3,nconstraints=1)
        gsm.opt_iters()
        if ORCA:
            os.system('rm temporcarun/*')

    #gsm.write_node_xyz('nodes_xyz_file1')