from cif2crystal import Cif2Crystal
from runcrystal import RunCrystal, RunProperties, NewRunProperties
from runqwalk import Crystal2QWalk, NewCrystal2QWalk,\
                     QWalkVarianceOptimize, QWalkRunVMC
from copy import deepcopy
import job_control as jc
import os
import json

import veritas as ver

results = []
baseroot = "agtest_"

default_job = jc.default_job_record("si/si.cif")
default_job['dft']['functional']['hybrid'] = 0
default_job['dft']['basis']=[0.2,2,2]
default_job['dft']['tolinteg']=[6,6,6,6,14]
default_job['dft']['kmesh']=[4,4,4]
default_job['dft']['spin_polarized']=False
default_job['qmc']['vmc']['optimizer'] = None
default_job['qmc']['nblocks'] = 10
default_job['qmc']['kpoints']='complex'

##############################################
# Simple Si calculation.
# Set defaults.
base = baseroot + "si_"
element_list = [
    Cif2Crystal(),
    RunCrystal(
      submitter=ver.LocalVeritasCrystalSubmitter(
        nn=1,np=8,time="20:00:00",queue="batch")),
    RunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch")),
    Crystal2QWalk(),
    #QWalkVarianceOptimize(
    #  submitter=ver.LocalVeritasQwalkSubmitter(
    #    nn=1,np=8,time="20:00:00",queue="batch")),
    QWalkRunVMC(
      submitter=ver.LocalVeritasQwalkSubmitter(
        nn=1,np=8,time="20:00:00",queue="batch"))
  ]

ref_job = deepcopy(default_job)
ref_job['control']['id'] = base+"ref"
results.append(jc.execute(ref_job,element_list))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(default_job)
test_job['control']['id'] = base+"test"
results.append(jc.execute(test_job,element_list))

##############################################
# Supercell Si calculation.
base = baseroot + "sup_"

element_list = [
    Cif2Crystal(),
    RunCrystal(
      submitter=ver.LocalVeritasCrystalSubmitter(
        nn=1,np=8,time="20:00:00",queue="batch")),
    RunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch")),
    Crystal2QWalk(),
    #QWalkVarianceOptimize(
    #  submitter=ver.LocalVeritasQwalkSubmitter(
    #    nn=1,np=8,time="20:00:00",queue="batch")),
    QWalkRunVMC(
      submitter=ver.LocalVeritasQwalkSubmitter(
        nn=1,np=8,time="20:00:00",queue="batch"))
  ]

default_job['supercell'] = [[2,0,0],[0,2,0],[0,0,2]]
default_job['dft']['kmesh']=[3,3,3] # ensures one complex.

ref_job = deepcopy(default_job)
ref_job['control']['id'] = base+"ref"
#results.append(jc.execute(ref_job,element_list))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(default_job)
test_job['control']['id'] = base+"test"
#results.append(jc.execute(test_job,element_list))

##############################################
# GaAs calculation (two different atoms).
base = baseroot + "gaas_"

element_list = [
    Cif2Crystal(),
    RunCrystal(
      submitter=ver.LocalVeritasCrystalSubmitter(
        nn=1,np=8,time="20:00:00",queue="batch")),
    RunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch")),
    Crystal2QWalk(),
    #QWalkVarianceOptimize(
    #  submitter=ver.LocalVeritasQwalkSubmitter(
    #    nn=1,np=8,time="20:00:00",queue="batch")),
    QWalkRunVMC(
      submitter=ver.LocalVeritasQwalkSubmitter(
        nn=1,np=8,time="20:00:00",queue="batch"))
  ]

default_job['cif'] = open("gaas/GaAs.cif",'r').read()
default_job['dft']['kmesh']=[3,3,3]

ref_job = deepcopy(default_job)
ref_job['control']['id'] = base+"ref"
results.append(jc.execute(ref_job,element_list))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(default_job)
test_job['control']['id'] = base+"test"
results.append(jc.execute(test_job,element_list))

###############################################
# Gather and export.
json.dump(results,open("test_results.json",'w'))
