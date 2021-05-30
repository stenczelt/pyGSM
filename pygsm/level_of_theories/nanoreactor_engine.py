# local application imports
from pygsm import utilities
from .base_lot import Lot

'''
'''

class nanoreactor_engine(Lot):

    def __init__(self,options):
        super(nanoreactor_engine,self).__init__(options)

        # can we do a check here?
        self.engine=options['job_data']['engine']
        self.nscffail = 0

    def run(self,geom,mult,ad_idx,runtype='gradient'):
        self.Gradients={}
        self.Energies = {}

        xyz = utilities.utilities.manage_xyz.xyz_to_np(geom) * utilities.utilities.units.ANGSTROM_TO_AU

        # Call the engine
        try:
            energy,gradient = self.engine.compute_gradient(xyz)
        except:
            # The calculation failed
            # set energy to a large number so the optimizer attempts to slow down
            print(" SCF FAILURE")
            self.nscffail+=1
            energy,gradient = 999, 0

            if self.nscffail>25:
                raise RuntimeError

        # Store the values in memory
        self._Energies[(mult,ad_idx)] = self.Energy(energy,'Hartree')
        self._Gradients[(mult,ad_idx)] = self.Gradient(gradient,'Hartree/Bohr')


if __name__=="__main__":
    from nanoreactor.engine import get_engine
    from nanoreactor.parsing import load_settings

    # read settings from name
    db, setting_name, settings = load_settings('refine.yaml', host='fire-05-30')

    # Create the nanoreactor TCPB engine
    engine_type=settings['engine'].pop('type')
    engine = get_engine(r.mol, engine_type=engine_type, **settings['engine'])


    # read in a geometry
    geom = utilities.manage_xyz.read_xyz('../../data/ethylene.xyz')
    xyz = utilities.manage_xyz.xyz_to_np(geom)

    # create the pygsm level of theory object
    test_lot = nanoreactor_engine(geom,job_data = {'engine',test_engine})

    # Test
    print("getting energy")
    print(test_lot.get_energy(xyz,1,0))

    print("getting grad")
    print(test_lot.get_gradient(xyz,1,0))



