from model_base import BaseModel
from darts.engines import sim_params, value_vector
from physics.physics_comp_sup import SuperPhysics
import numpy as np
from physics.properties_basic import *
from physics.property_container import *

# Model class creation here!
class Model(BaseModel):
    def __init__(self):
        # Call base class constructor
        super().__init__()

        self.zero = 1e-12
        """Physical properties"""
        # Create property containers:
        components_name = ['CO2', 'C1', 'H2O']
        Mw = [44.01, 16.04, 18.015]
        self.property_container = property_container(phases_name=['gas', 'wat'],
                                                     components_name=components_name,
                                                     Mw=Mw, min_z=self.zero / 10)
        self.components = self.property_container.components_name
        self.phases = self.property_container.phases_name

        """ properties correlations """
        self.property_container.flash_ev = Flash(self.components, [4, 2, 1e-1], self.zero)
        self.property_container.density_ev = dict([('gas', Density(compr=1e-3, dens0=200)),
                                                   ('wat', Density(compr=1e-5, dens0=600))])
        self.property_container.viscosity_ev = dict([('gas', ViscosityConst(0.05)),
                                                     ('wat', ViscosityConst(0.5))])
        self.property_container.rel_perm_ev = dict([('gas', PhaseRelPerm("gas")),
                                                    ('wat', PhaseRelPerm("oil"))])

        """ Activate physics """
        self.physics = SuperPhysics(self.property_container, self.timer, n_points=200, min_p=1, max_p=400,
                                    min_z=self.zero/10, max_z=1-self.zero/10)

        self.inj_comp = [1.0 - 2 * self.zero, self.zero]
        self.ini_comp = [0.01, 0.2]
        self.pressure_ini = 100

        # Some newton parameters for non-linear solution:
        self.params.first_ts = 1e-5
        self.params.mult_ts = 2
        self.params.max_ts = 5

        self.params.tolerance_newton = 1e-2
        self.params.tolerance_linear = 1e-3
        self.params.max_i_newton = 20
        self.params.max_i_linear = 30
        self.params.newton_type = sim_params.newton_local_chop

    def set_boundary_conditions(self):
        for i, w in enumerate(self.reservoir.wells):
            if 'INJ' in w.name:
                w.control = self.physics.new_rate_inj(0, self.inj_comp, 0)
            else:
                w.control = self.physics.new_rate_prod(0, 1)

    def set_wells(self):
        for i, w in enumerate(self.reservoir.wells):
            if w.name == 'INJ017':
                w.control = self.physics.new_bhp_inj(120, self.inj_comp)
            elif w.name == 'PROD024A' or w.name == 'PROD009':
                w.control = self.physics.new_bhp_prod(80)
            else:
                w.control = self.physics.new_rate_prod(0, 1)

    def set_op_list(self):
        self.op_num = np.array(self.reservoir.mesh.op_num, copy=False)
        n_res = self.reservoir.mesh.n_res_blocks
        self.op_num[n_res:] = 1
        self.op_list = [self.physics.acc_flux_itor, self.physics.acc_flux_w_itor]

    def export_pro_vtk(self, file_name='co2_out'):
        Xn = np.array(self.physics.engine.X, copy=False)
        P = Xn[0:self.reservoir.nb * 3:3]
        z1 = Xn[1:self.reservoir.nb * 3:3]
        z2 = Xn[2:self.reservoir.nb * 3:3]

        sg = np.zeros(len(P))
        sw = np.zeros(len(P))

        for i in range(len(P)):
            values = value_vector([0] * self.physics.n_ops)
            state = value_vector((P[i], z1[i], z2[i]))
            self.physics.property_itor.evaluate(state, values)
            sg[i] = values[0]
            sw[i] = 1 - sg[i]

        self.export_vtk(file_name, local_cell_data={'GasSat': sg, 'WatSat': sw})

