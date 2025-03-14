'''
pull data and run curve fitting code
'''

import pandas as pd
import synapseclient as sc
import os
import subprocess
syn = sc.login()
##this is a pain, need to move to file
curtabs={'NF0019_T1':"syn64954302",\
        'NF0017_T3':"syn64954273",\
        'NF0023_T2':"syn64938448",\
        'NF0022_T1':"syn64938446",\
        'NF0020_T1':"syn64938445",\
        'NF0019_T1':"syn64938444",\
        'NF0018_T1':"syn64938443"}

  
###pull files
patlist=[]
for pat, id in curtabs.items():
  tab = syn.tableQuery('select * from '+id).asDataFrame()
  tab["improve_sample_id"] = pat
  tab = tab.reset_index()
  print(tab)
  patlist.append(tab)
  
fulltab = pd.concat(patlist)
fulltab['DOSE']=fulltab.Concentration_uM+0.0001
fulltab['GROWTH']=fulltab.Viability_percentage
fulltab['time']=120
fulltab['time_unit']='hours'
fulltab['study']='cnfPDO'
fulltab['source']='synapse'
##mutate the values create new columns
ncols=['DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit']
fulltab = fulltab[ncols]
##change file headers to DOSE/RESPONSE values needed by other script
fulltab.to_csv('cohort1_drug_response.tsv',sep='\t')

###store on synapse
syn.store(sc.File('cohort_1_drug_response.tsv',parentId='syn65471813'))

script='https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/build/utils/fit_curve.py'
subprocess.run(['wget',script])
subprocess.run(['python','fit_curve.py','--input','drug_response.tsv','--output','cnfDrugOutput'])
##rename file
os.system('mv cnfDrugOutput.0 cohort1_curves.tsv')
syn.store(sc.File('cohort1_curves.tsv',parentId='syn65471813'))


#store on synapse
