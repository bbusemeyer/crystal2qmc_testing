from __future__ import print_function
from cif2crystal import Cif2Crystal
from xyz2crystal import Xyz2Crystal
from runcrystal import RunCrystal, RunProperties, NewRunProperties
from runqwalk import Crystal2QWalk, NewCrystal2QWalk,\
                     QWalkVarianceOptimize, QWalkRunVMC
from copy import deepcopy
import job_control as jc
import os
import json
import check_results

import veritas as ver

results = []
baseroot = "agtest_"

default_job = jc.default_job_record("si/si.cif")
default_job['dft']['functional']['hybrid'] = 0
default_job['dft']['basis']=[0.2,2,2]
default_job['dft']['tolinteg']=[6,6,6,6,14]
default_job['dft']['edifftol'] = 6
default_job['dft']['kmesh']=[4,4,4]
default_job['dft']['spin_polarized']=False
default_job['qmc']['vmc']['optimizer'] = None
default_job['qmc']['vmc']['nblock'] = 100
default_job['qmc']['vmc']['target_error'] = 0.01
default_job['qmc']['kpoints']='complex'

##############################################
# Simple Si calculation.
# Set defaults.
checking_this = True
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
if checking_this:
  results.append(jc.execute(ref_job,element_list))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(default_job)
test_job['control']['id'] = base+"test"
test_job['dft']['restart_from'] = "../"+base+"ref"+"/fort.79"
if checking_this:
  results.append(jc.execute(test_job,element_list))

##############################################
# Supercell Si calculation.
checking_this = False
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

cur_job = deepcopy(default_job)
cur_job['supercell'] = [[2,0,0],[0,2,0],[0,0,2]]
cur_job['dft']['kmesh']=[3,3,3] # ensures one complex.

ref_job = deepcopy(cur_job)
ref_job['control']['id'] = base+"ref"
if checking_this:
  results.append(jc.execute(ref_job,element_list[:3]))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(cur_job)
test_job['control']['id'] = base+"test"
test_job['dft']['restart_from'] = "../"+base+"ref"+"/fort.79"
if checking_this:
  results.append(jc.execute(test_job,element_list[:3]))

##############################################
# GaAs calculation (two different atoms).
checking_this = True
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

cur_job = deepcopy(default_job)
cur_job['cif'] = open("gaas/GaAs.cif",'r').read()
cur_job['dft']['kmesh']=[3,3,3]

ref_job = deepcopy(cur_job)
ref_job['control']['id'] = base+"ref"
if checking_this:
  results.append(jc.execute(ref_job,element_list))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(cur_job)
test_job['control']['id'] = base+"test"
test_job['dft']['restart_from'] = "../"+base+"ref"+"/fort.79"
if checking_this:
  results.append(jc.execute(test_job,element_list))

##############################################
# Fe calculation (Spin enabled).
checking_this = True
base = baseroot + "fe_"

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

cur_job = deepcopy(default_job)
cur_job['cif'] = open("Fe.cif",'r').read()
cur_job['dft']['kmesh']=[3,3,3]
cur_job['dft']['spin_polarized']=True
cur_job['dft']['initial_spin'] = [1]
cur_job['qmc']['vmc']['target_error'] = 0.02
cur_job['total_spin'] = 2

ref_job = deepcopy(cur_job)
ref_job['control']['id'] = base+"ref"
if checking_this:
  results.append(jc.execute(ref_job,element_list))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(cur_job)
test_job['control']['id'] = base+"test"
test_job['dft']['restart_from'] = "../"+base+"ref"+"/fort.79"
if checking_this:
  results.append(jc.execute(test_job,element_list))

##############################################
# KMnF3 calculation (Spin enabled).
checking_this = True
base = baseroot + "kmnf3_"

element_list = [
    Cif2Crystal(),
    RunCrystal(
      submitter=ver.LocalVeritasCrystalSubmitter(
        nn=1,np=8,time="120:00:00",queue="batch")),
    RunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="120:00:00",queue="batch")),
    Crystal2QWalk(),
    #QWalkVarianceOptimize(
    #  submitter=ver.LocalVeritasQwalkSubmitter(
    #    nn=1,np=8,time="20:00:00",queue="batch")),
    QWalkRunVMC(
      submitter=ver.LocalVeritasQwalkSubmitter(
        nn=1,np=8,time="120:00:00",queue="batch"))
  ]

cur_job = deepcopy(default_job)
cur_job['cif'] = open("KMnF3.cif",'r').read()
cur_job['dft']['kmesh']=[3,3,3]
cur_job['dft']['spin_polarized']=True
cur_job['dft']['initial_spin'] = [0,1,0,0,0]
cur_job['qmc']['vmc']['target_error'] = 0.05
cur_job['total_spin'] = 5

ref_job = deepcopy(cur_job)
ref_job['control']['id'] = base+"ref"
if checking_this:
  results.append(jc.execute(ref_job,element_list))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="120:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(cur_job)
test_job['control']['id'] = base+"test"
test_job['dft']['restart_from'] = "../"+base+"ref"+"/fort.79"
if checking_this:
  results.append(jc.execute(test_job,element_list))

##############################################
# MnO calculation (Antiferromagnetic, expensive).
checking_this = False
base = baseroot + "mno_"

element_list = [
    Cif2Crystal(),
    RunCrystal(
      submitter=ver.LocalVeritasCrystalSubmitter(
        nn=1,np=8,time="120:00:00",queue="batch")),
    RunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="120:00:00",queue="batch")),
    Crystal2QWalk(),
    #QWalkVarianceOptimize(
    #  submitter=ver.LocalVeritasQwalkSubmitter(
    #    nn=1,np=8,time="20:00:00",queue="batch")),
    QWalkRunVMC(
      submitter=ver.LocalVeritasQwalkSubmitter(
        nn=1,np=8,time="120:00:00",queue="batch"))
  ]

cur_job = deepcopy(default_job)
cur_job['cif'] = open("mno/MnO.cif",'r').read()
cur_job['supercell'] = [[2,0,0],[0,2,0],[0,0,2]]
cur_job['dft']['initial_spin'] = [1,1,1,-1,1,-1,-1,-1,0,0,0,0,0,0,0,0]
cur_job['dft']['kmesh']=[3,3,3]
cur_job['dft']['spin_polarized']=True
cur_job['qmc']['vmc']['target_error'] = 0.05

ref_job = deepcopy(cur_job)
ref_job['control']['id'] = base+"ref"
if checking_this:
  results.append(jc.execute(ref_job,element_list))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="120:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(cur_job)
test_job['control']['id'] = base+"test"
test_job['dft']['restart_from'] = "../"+base+"ref"+"/fort.79"
if checking_this:
  results.append(jc.execute(test_job,element_list))


##############################################
# CO (molecule).
checking_this = True
base = baseroot + "co_"

element_list = [
    Xyz2Crystal(),
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

del cur_job['cif']
cur_job['xyz'] = open("CO.xyz",'r').read()
cur_job['qmc']['vmc']['target_error'] = 0.005

ref_job = deepcopy(cur_job)
ref_job['control']['id'] = base+"ref"
if checking_this:
  results.append(jc.execute(ref_job,element_list))

element_list[2] = \
    NewRunProperties(
      submitter=ver.LocalVeritasPropertiesSubmitter(
        nn=1,np=1,time="20:00:00",queue="batch"))
element_list[3] = NewCrystal2QWalk()

test_job = deepcopy(cur_job)
test_job['control']['id'] = base+"test"
test_job['dft']['restart_from'] = "../"+base+"ref"+"/fort.79"
if checking_this:
  results.append(jc.execute(test_job,element_list))

##############################################
# Graphene (2-d system)
checking_this = True
base = baseroot + "graphene_"
element_list = [
    QWalkRunVMC(
      submitter=ver.LocalVeritasQwalkSubmitter(
        nn=1,np=8,time="20:00:00",queue="batch"))
  ]

cur_job['qmc']['vmc']['target_error'] = 0.01

ref_job = deepcopy(cur_job)
ref_job['control']['id'] = base+"ref"
if checking_this:
  results.append(jc.execute(ref_job,element_list))

test_job = deepcopy(cur_job)
test_job['control']['id'] = base+"test"
test_job['dft']['restart_from'] = "../"+base+"ref"+"/fort.79"
if checking_this:
  results.append(jc.execute(test_job,element_list))

###############################################
# Gather and export.
json.dump(results,open("test_results.json",'w'))

print("~~~~~~~~~~~~~~~~~~~~~~")
print("Checking results...")
try:
  check_results.perform_check("test_results.json")
except KeyError:
  print("KeyError occured, maybe no calculations have finished yet?")
except IndexError:
  print("IndexError occured, maybe only one type of calculation has finished yet?")
print("======================")
