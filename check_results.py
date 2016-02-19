from __future__ import print_function
import pandas as pd
import data_processing as dp

def process_ids(idstr):
    return pd.Series(idstr.split('_')[1:],['material','type'])

def compare_results(df,tol=3.0):
  count = 0
  allcheck = True
  kcheck   = {}
  dres  = df.loc[df['type']=='ref', 'results'].iloc[0]
  dtest = df.loc[df['type']=='test','results'].iloc[0]
  for key in dres.keys():
    count += 1
    check = dtest[key][0] - dres[key][0] \
        < tol*(dres[key][1]**2 + dtest[key][1]**2)**.5
    kcheck[key] = check
    allcheck   &= check
    if not check:
      print("Fail %s k%d: %f > %f"%(
        df['material'].iloc[0],
        key,
        dtest[key][0] - dres[key][0],
        tol*(dres[key][1]**2 + dtest[key][1]**2)**.5))
  #print "Total count %d."%count
  refegy  = df.loc[df['type']=='ref', 'vmc_energy'].iloc[0]
  testegy = df.loc[df['type']=='test', 'vmc_energy'].iloc[0]
  referr  = df.loc[df['type']=='ref', 'vmc_error'].iloc[0]
  testerr = df.loc[df['type']=='test', 'vmc_error'].iloc[0]
  avgcheck = testegy-refegy < tol*(testerr**2 + referr**2)**.5
  return pd.Series([kcheck,allcheck,avgcheck],['kcheck','allcheck','avgcheck'])

def perform_check(inpjson="test_results.json"):
  rawdf,alldf = dp.format_autogen(inpjson)
  alldf = alldf.join(alldf['id'].apply(process_ids))
  alldf = dp.kavergage_qmc(alldf,qmc_type='vmc')
  alldf['results'] = alldf.loc[alldf['results'].notnull(),'results']\
      .apply(dp.format_results)
  checkdf = alldf[alldf['results'].notnull()]\
      .groupby('material')\
      .apply(compare_results)
  if checkdf['allcheck'].all() and checkdf['avgcheck'].all():
    print("Clear pass!")
  else:
    print("Potential problem.")
    if not checkdf['allcheck'].all():
      print("Some k-points don't match:")
      print(checkdf.loc[~alldf['allcheck'],:])
    if not checkdf['avgcheck'].all():
      print("Average doesn't match:")
      print(checkdf.loc[~checkdf['avgcheck'],:])
  if alldf['results'].isnull().any():
    print("Although, some are not done.")

if __name__ == "__main__":
  perform_check()
