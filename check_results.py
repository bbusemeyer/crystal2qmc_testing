from __future__ import print_function
import pandas as pd
import data_processing as dp

def process_ids(idstr):
    return pd.Series(idstr.split('_')[1:],['material','type'])

def compare_results(df,tol=2.5):
    count = 0
    allcheck = True
    kcheck   = {}
    dres  = df.loc[df['type']=='ref', 'results'].iloc[0]
    dtest = df.loc[df['type']=='test','results'].iloc[0]
    for key in dres.keys():
        count += 1
        check = dtest[key][0] - dres[key][0] < tol*(dres[key][1]**2 + dtest[key][1]**2)**.5
        kcheck[key] = check
        allcheck   &= check
        #if not check:
        #    print "Fail %d: %f > %f"%(key,dtest[key][0] - dres[key][0],tol*(dres[key][1]**2 + dtest[key][1]**2)**.5)
    #print "Total count %d."%count
    return pd.Series([kcheck,allcheck],['kcheck','allcheck'])

def perform_check(inpjson="test_results.json"):
  rawdf,alldf = dp.format_autogen(inpjson)
  alldf = alldf.join(alldf['id'].apply(process_ids))
  alldf['results'] = alldf.loc[alldf['results'].notnull(),'results']\
      .apply(dp.format_results)
  checkdf = alldf[alldf['results'].notnull()]\
      .groupby('material')\
      .apply(compare_results)
  if checkdf['allcheck'].all():
    print("Clear pass!")
  else:
    print("Potential problem.")
  if alldf['results'].isnull().any():
    print("Although, some are not done.")

if __name__ == "__main__":
  perform_check()
