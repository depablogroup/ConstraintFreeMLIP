from ase import units
from ase.io import read, write
from ase.calculators.cp2k import CP2K


# Temperature and timestep
T = 300
dt = 1 * units.fs


inp = """
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME  BASIS_SET
    POTENTIAL_FILE_NAME  GTH_POTENTIALS
    &MGRID
      CUTOFF 550
      REL_CUTOFF 50 # was 60
    &END MGRID
    &QS
      METHOD GPW
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
      EPS_DEFAULT 1.0E-10
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 500
      &OT
        PRECONDITIONER FULL_SINGLE_INVERSE
        MINIMIZER CG
        LINESEARCH 3PNT
      &END
      &OUTER_SCF
        EPS_SCF 1.0E-5
        MAX_SCF 200
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        DISPERSION_FUNCTIONAL PAIR_POTENTIAL
        &PAIR_POTENTIAL
         TYPE DFTD3
         PARAMETER_FILE_NAME dftd3.dat
         REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &KIND C
      ELEMENT   C
      BASIS_SET DZVP-GTH-PBE
      POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND H
      ELEMENT H
      BASIS_SET DZVP-GTH-PBE
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND N
      ELEMENT N
      BASIS_SET DZVP-GTH-PBE
      POTENTIAL GTH-PBE-q5
    &END KIND
    &KIND O
      ELEMENT O
      BASIS_SET DZVP-GTH-PBE
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
"""

loops=open('potenersese.dat','a')
dimer = read('config.extxyz', index="1")
dimer.set_cell([32.168, 32.168, 32.168])
dimer.set_pbc([1,1,1])

calc = CP2K(basis_set=None,
        basis_set_file=None,
        max_scf=None,
        cutoff=None,
        force_eval_method=None,
        potential_file=None,
        poisson_solver=None,
        pseudo_potential=None,
        stress_tensor=False,
        xc=None,
        inp=inp)
dimer.set_calculator(calc)
e = dimer.get_potential_energy()
write('train.extxyz', dimer,append=True)
loops.write(str(e)+'\n')
loops.close()

