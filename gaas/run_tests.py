import cif2crystal
import runcrystal
import runqwalk
import copy
import job_control as jc
import os
import json

import veritas

element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
element_list.append(runcrystal.RunCrystal(
  submitter=veritas.LocalVeritasCrystalSubmitter(
    nn=1,np=8,time="0:20:00",queue="batch")))
element_list.append(runcrystal.RunProperties(
  submitter=veritas.LocalVeritasPropertiesSubmitter(
    nn=1,np=1,time="1:00:00",queue="batch")))
element_list.append(runqwalk.Crystal2QWalk())
element_list.append(runqwalk.QWalkVarianceOptimize(
  submitter=veritas.LocalVeritasQwalkSubmitter(
    nn=1,np=8,time="0:10:00",queue="batch")))
element_list.append(runqwalk.QWalkEnergyOptimize(
  submitter=veritas.LocalVeritasQwalkSubmitter(
    nn=1,np=8,time="0:10:00",queue="batch")))
element_list.append(runqwalk.QWalkRunDMC(
  submitter=veritas.LocalVeritasBundleQwalkSubmitter(
    nn=1,np=8,time="0:10:00",queue="batch")))

default_job=jc.default_job_record("GaAs.cif")
default_job['dft']['kmesh'] = [2,2,2]
default_job['dft']['functional']['hybrid'] = 0
default_job['dft']['tolinteg'] = [6,6,6,6,12]
default_job['dft']['basis']=[0.2,2,2]
default_job['dft']['maxcycle'] = 100
default_job['dft']['fmixing'] = 80
default_job['dft']['edifftol'] = 6
default_job['dft']['broyden'] = [0.1,60,20]
default_job['qmc']['variance_optimize']['reltol']=0.1
default_job['qmc']['variance_optimize']['abstol']=10
default_job['qmc']['dmc']['save_trace'] = False
default_job['qmc']['dmc']['nblock']=5
default_job['qmc']['dmc']['target_error']=0.1
default_job['total_spin'] = 0

name = "referece"
job_record = copy.deepcopy(default_job)
job_record['control']['id']=name
jc.execute(job_record,element_list)
