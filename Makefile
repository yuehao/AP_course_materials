html_maker = jupyter nbconvert
flag = --to html_with_toclenvs

topics= magnets chromatic_effect electron_storage_ring fringe_effect linear_transverse_coupling linear_transverse_motion longitudinal_dynamics off_momentum_orbit perturbed_linear_transverse_motion symplecticity synchrotron_radiation twiss_param undulator_radiation DBA FODO


.PHONY: all $(topics)

all: $(topics)

magnets: magnets/magnets.ipynb
	$(html_maker) $(flag) $^

chromatic_effect: chromatic_effect/chromatic_effect.ipynb
	$(html_maker) $(flag) $^

electron_storage_ring:  electron_storage_ring/electron_storage_ring.ipynb
	$(html_maker) $(flag) $^

fringe_effect: fringe_effect/fringe_effect.ipynb
	$(html_maker) $(flag) $^

linear_transverse_coupling:linear_transverse_coupling/linear_transverse_coupling.ipynb
	$(html_maker) $(flag) $^

linear_transverse_motion:linear_transverse_motion/linear_transverse_motion.ipynb
	$(html_maker) $(flag) $^

longitudinal_dynamics: longitudinal_dynamics/longitudinal_dynamics.ipynb
	$(html_maker) $(flag) $^

off_momentum_orbit: off_momentum_orbit/off_momentum_orbit.ipynb
	$(html_maker) $(flag) $^

perturbed_linear_transverse_motion:perturbed_linear_transverse_motion/perturbed_linear_transverse_motion.ipynb
	$(html_maker) $(flag) $^

symplecticity:symplecticity/symplecticity.ipynb
	$(html_maker) $(flag) $^

synchrotron_radiation:synchrotron_radiation/synchrotron_radiation.ipynb
	$(html_maker) $(flag) $^

twiss_param: twiss_param/twiss_param.ipynb
	$(html_maker) $(flag) $^

undulator_radiation: undulator_radiation/undulator_radiation.ipynb
	$(html_maker) $(flag) $^

DBA: examples/DBA/DBA.ipynb
	$(html_maker) $(flag) $^

FODO: examples/FODO/FODO_cell.ipynb
	$(html_maker) $(flag) $^
