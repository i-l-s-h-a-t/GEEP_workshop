from physics.geothermal import Geothermal
from darts.models.physics.iapws.iapws_property import *
from darts.models.physics.iapws.custom_rock_property import *
from model_base import BaseModel
from darts.engines import value_vector
import numpy as np
from physics.property_container import *


class Model(BaseModel):

    def __init__(self, n_points=100):
        # call base class constructor
        super().__init__()

        self.hcap = np.array(self.reservoir.mesh.heat_capacity, copy=False)
        self.conduction = np.array(self.reservoir.mesh.rock_cond, copy=False)
        self.hcap.fill(2200)
        self.conduction.fill(181.44)

        # Create property containers:
        self.property_container = model_properties(phases_name=['water', 'steam', 'temperature', 'energy'],
                                                   components_name=['H2O'])

        # Define properties in property_container (IAPWS is the default property package for Geothermal in DARTS)
        # Users can define their custom properties in custom_properties.py; several property examples are defined there.
        self.rock = [value_vector([1, 0, 273.15])]
        self.property_container.temp_ev = iapws_temperature_evaluator()
        self.property_container.enthalpy_ev = dict([('water', iapws_water_enthalpy_evaluator()),
                                                    ('steam', iapws_steam_enthalpy_evaluator())])
        self.property_container.saturation_ev = dict([('water', iapws_water_saturation_evaluator()),
                                                    ('steam', iapws_steam_saturation_evaluator())])
        self.property_container.rel_perm_ev = dict([('water', iapws_water_relperm_evaluator()),
                                                    ('steam', iapws_steam_relperm_evaluator())])
        self.property_container.density_ev = dict([('water', iapws_water_density_evaluator()),
                                                   ('steam', iapws_steam_density_evaluator())])
        self.property_container.viscosity_ev = dict([('water', iapws_water_viscosity_evaluator()),
                                                     ('steam', iapws_steam_viscosity_evaluator())])
        self.property_container.saturation_ev = dict([('water', iapws_water_saturation_evaluator()),
                                                      ('steam', iapws_steam_saturation_evaluator())])

        self.property_container.rock_compaction_ev = custom_rock_compaction_evaluator(self.rock)
        self.property_container.rock_energy_ev = custom_rock_energy_evaluator(self.rock)

        self.physics = Geothermal(property_container=self.property_container, timer=self.timer, n_points=1000, min_p=1,
                                  max_p=500, min_e=100, max_e=55000, grav=False)

        self.params.first_ts = 1e-5
        self.params.mult_ts = 2
        self.params.max_ts = 100

        # Newton tolerance is relatively high because of L2-norm for residual and well segments
        self.params.tolerance_newton = 1e-3
        self.params.tolerance_linear = 1e-3
        self.params.max_i_newton = 20
        self.params.max_i_linear = 30
        self.runtime = 3650

    def set_initial_conditions(self):
        mesh = self.reservoir.mesh

        pres_init = 200
        temp_init = 273 + 60
        pressure = np.array(mesh.pressure, copy=False)
        pressure[:] = pres_init

        enthalpy = np.array(mesh.enthalpy, copy=False)

        state = value_vector([pres_init, 0])
        E = iapws_total_enthalpy_evalutor(temp_init)
        enthalpy_init = E.evaluate(state)
        enthalpy[:] = enthalpy_init


    def set_boundary_conditions(self):
        for i, w in enumerate(self.reservoir.wells):
            if "INJ" in w.name:
                w.control = self.physics.new_rate_water_inj(0, 308)
            else:
                w.control = self.physics.new_rate_water_prod(0)

    def set_wells(self):
        for i, w in enumerate(self.reservoir.wells):
            if "INJ" in w.name:
                w.control = self.physics.new_bhp_water_inj(220, 308)
            else:
                w.control = self.physics.new_bhp_prod(180)

    def enthalpy_to_temperature(self, data):
        from darts.models.physics.iapws.iapws_property_vec import _Backward1_T_Ph_vec
        data_len = int(len(data) / 2)
        T = np.zeros(data_len)
        T[:] = _Backward1_T_Ph_vec(data[::2] / 10, data[1::2] / 18.015)
        return T

    def export_pro_vtk(self, file_name='geothermal_out'):
        from darts.models.physics.iapws.iapws_property_vec import _Backward1_T_Ph_vec
        nb = self.reservoir.nb
        T = np.zeros(nb)
        X = np.array(self.physics.engine.X)
        T[:] = _Backward1_T_Ph_vec(X[:nb * 2:2] / 10, X[1:nb * 2:2] / 18.015)
        self.export_vtk(file_name=file_name, local_cell_data={'temp': T})

class model_properties(property_container):
    def __init__(self, phases_name, components_name):
        # Call base class constructor
        self.nph = len(phases_name)
        Mw = np.ones(self.nph)
        super().__init__(phases_name, components_name, Mw)

        # remove the virtual phase from the parent class
        self.dens = np.zeros(self.nph-2)
        self.sat = np.zeros(self.nph-2)
        self.mu = np.zeros(self.nph-2)
        self.kr = np.zeros(self.nph-2)
        self.enthalpy = np.zeros(self.nph-2)

    def evaluate(self, state):
        vec_state_as_np = np.asarray(state)

        for j in range(self.nph-2):
            self.enthalpy[j] = self.enthalpy_ev[self.phases_name[j]].evaluate(state)
            self.dens[j] = self.density_ev[self.phases_name[j]].evaluate(state)
            self.sat[j] = self.saturation_ev[self.phases_name[j]].evaluate(state)
            self.kr[j] = self.rel_perm_ev[self.phases_name[j]].evaluate(state)
            self.mu[j] = self.viscosity_ev[self.phases_name[j]].evaluate(state)

        self.temp = self.temp_ev.evaluate(state)
        self.rock_compaction = self.rock_compaction_ev.evaluate(state)
        self.rock_int_energy = self.rock_energy_ev.evaluate(state)

        return self.enthalpy, self.dens, self.sat, self.kr, self.mu, self.temp, self.rock_compaction, self.rock_int_energy



