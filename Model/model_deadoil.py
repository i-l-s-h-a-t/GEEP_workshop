from physics.physics_comp_sup import SuperPhysics
from darts.engines import value_vector
from model_base import BaseModel
import numpy as np
from physics.property_container import *
from physics.properties_dead_oil import *


class Model(BaseModel):

    def __init__(self):
        # call base class constructor
        super().__init__()

        """Physical properties"""
        self.pvt = 'physics_do.in'
        self.zero = 1e-13
        self.property_container = model_properties(phases_name=['water', 'oil'], components_name=['w', 'o'],
                                                   pvt=self.pvt, min_z=self.zero / 10)

        # Define property evaluators based on custom properties
        self.flash_ev = []
        self.property_container.density_ev = dict([('water', DensityWat(self.pvt)),
                                                   ('oil', DensityOil(self.pvt))])
        self.property_container.viscosity_ev = dict([('water', ViscoWat(self.pvt)),
                                                     ('oil', ViscoOil(self.pvt))])
        self.property_container.rel_perm_ev = dict([('water', WatRelPerm(self.pvt)),
                                                    ('oil', OilRelPerm(self.pvt))])
        self.property_container.capillary_pressure_ev = CapillarypressurePcow(self.pvt)

        self.property_container.rock_compress_ev = RockCompactionEvaluator(self.pvt)

        # create physics
        self.thermal = 0
        self.physics = SuperPhysics(self.property_container, self.timer, n_points=400, min_p=0, max_p=1000,
                                    min_z=self.zero, max_z=1 - self.zero, thermal=self.thermal)
        self.params.first_ts = 1e-3
        self.params.mult_ts = 2
        self.params.max_ts = 100
        self.params.tolerance_newton = 1e-3
        self.params.tolerance_linear = 1e-4

        self.runtime = 300
        self.ini_comp = value_vector([0.2357])
        self.inj_comp = value_vector([1 - self.zero])
        self.pressure_ini = 400

        self.timer.node["initialization"].stop()

    def print_array_stat(self, name, arr):
        print(name, 'MIN=', arr.min(), 'MEAN=', arr.mean(), 'MAX=', arr.max())


    def set_initial_conditions(self):
        """ Uniform Initial conditions for dead oil
            Set pressure by depth gradient (start from 1 bar at zero depth)
            Fill composition from initial water saturation self.sw
        """
        mesh = self.reservoir.mesh
        nb = mesh.n_blocks

        depth = self.depth[self.actnum>0]
        #depth = np.array(mesh.depth, copy=True)
        self.print_array_stat('depth', depth)

        # set initial pressure
        pressure_grad = 100  # bars per meter
        pressure = np.array(mesh.pressure, copy=False)
        pressure[:nb - 2 * len(self.reservoir.wells)] = depth / 1000. * pressure_grad + 1.
        self.print_array_stat('initial_pressure', pressure)

        # set initial composition
        nc = self.property_container.nc

        mesh.composition.resize(nb * (nc - 1))
        composition = np.array(mesh.composition, copy=False)
        # calculate composition from initial saturation
        # use surface densities. However, it is better to use densities at reservoir pressure
        sw = self.generate_sw(swl=0.2, woc_depth=3150, file_name='swat.inc')
        #sw = load_single_keyword(self.prop_filename, 'SWAT')
        self.print_array_stat('initial_Sw', sw)
        sw = sw[self.actnum > 0]
        so = 1 - sw
        dens_w = self.property_container.surf_wat_dens
        dens_o = self.property_container.surf_oil_dens
        composition[:nb - 2 * len(self.reservoir.wells)] = (sw * dens_w) / (sw * dens_w + so * dens_o) - self.zero
        self.print_array_stat('initial_composition', composition)

    def set_boundary_conditions(self):
        '''
        set well control by rate = 0
        '''
        for i, w in enumerate(self.reservoir.wells):
            if "INJ" in w.name:
                w.control = self.physics.new_rate_inj(0, self.inj_comp, 0)
            else:
                w.control = self.physics.new_rate_prod(0, 1)

    def set_wells(self):
        '''
        set well control by BHP
        '''
        for i, w in enumerate(self.reservoir.wells):
            if "INJ" in w.name:
                w.control = self.physics.new_rate_inj(2000, self.inj_comp, 0)
                w.constraint = self.physics.new_bhp_inj(550, self.inj_comp)
            else:
                w.control = self.physics.new_bhp_prod(350)

    def export_pro_vtk(self, file_name='deadoil_out'):
        Xn = np.array(self.physics.engine.X, copy=False)
        P = Xn[0:self.reservoir.nb * 2:2]
        z1 = Xn[1:self.reservoir.nb * 2:2]

        so = np.zeros(len(P))
        sw = np.zeros(len(P))

        for i in range(len(P)):
            values = value_vector([0] * self.physics.n_ops)
            state = value_vector((P[i], z1[i]))
            self.physics.property_itor.evaluate(state, values)
            sw[i] = values[0]
            so[i] = 1 - sw[i]

        self.export_vtk(file_name, local_cell_data={'OilSat': so, 'WatSat': sw})

    def generate_sw(self, swl, woc_depth, file_name='swat.inc'):
        '''
        fill water saturation cube with value=swl for cells with depth <= woc_depth,
                                   with value=1 otherwise
        :param swl: connate water saturation
        :param woc_depth: water oil contact depth
        :param file_name: file to write SWAT cube (GRDECL)
        :return: swat 1d-array
        '''
        nb = self.reservoir.nx * self.reservoir.ny * self.reservoir.nz
        sw = np.zeros(nb) + swl
        sw[self.depth > woc_depth] = 1.0
        save_few_keywords(file_name, ['SWAT'], [sw])
        return sw


class model_properties(property_container):
    def __init__(self, phases_name, components_name, pvt, min_z=1e-11):
        # Call base class constructor
        self.nph = len(phases_name)
        Mw = np.ones(self.nph)
        super().__init__(phases_name, components_name, Mw, min_z)
        self.x = np.zeros((self.nph, self.nc))
        self.pvt = pvt
        self.surf_dens = get_table_keyword(self.pvt, 'DENSITY')[0]
        self.surf_oil_dens = self.surf_dens[0]
        self.surf_wat_dens = self.surf_dens[1]

    def evaluate(self, state):
        """
        Class methods which evaluates the state operators for the element based physics
        :param state: state variables [pres, comp_0, ..., comp_N-1]
        :param values: values of the operators (used for storing the operator values)
        :return: updated value for operators, stored in values
        """
        # Composition vector and pressure from state:
        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]

        zc = np.append(vec_state_as_np[1:], 1 - np.sum(vec_state_as_np[1:]))

        self.clean_arrays()
        # two-phase flash - assume water phase is always present and water component last
        for i in range(self.nph):
            self.x[i, i] = 1

        ph = [0, 1]

        for j in ph:
            M = 0
            # molar weight of mixture
            for i in range(self.nc):
                M += self.Mw[i] * self.x[j][i]
            self.dens[j] = self.density_ev[self.phases_name[j]].evaluate(state)  # output in [kg/m3]
            self.dens_m[j] = self.dens[j] / M
            self.mu[j] = self.viscosity_ev[self.phases_name[j]].evaluate(state)  # output in [cp]

        self.nu = zc
        self.compute_saturation(ph)

        for j in ph:
            self.kr[j] = self.rel_perm_ev[self.phases_name[j]].evaluate(self.sat[0])
            self.pc[j] = 0

        return self.sat, self.x, self.dens, self.dens_m, self.mu, self.kr, self.pc, ph

    def evaluate_at_cond(self, pressure, zc):

        self.sat[:] = 0

        state = value_vector([1, 0])

        ph = [0, 1]
        for j in ph:
            self.dens_m[j] = self.density_ev[self.phases_name[j]].evaluate(state)

        self.dens_m = [self.surf_wat_dens, self.surf_oil_dens]  # to match DO based on PVT

        self.nu = zc
        self.compute_saturation(ph)

        return self.sat, self.dens_m
