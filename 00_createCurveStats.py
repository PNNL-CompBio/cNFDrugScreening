'''
pull data and run curve fitting code
'''

import pandas as pd
import synapseclient as sc
import os
import subprocess
syn = sc.login()
##this is a pain, need to move to file

filelist = syn.tableQuery("select id,individualID,specimenID from syn51301431 where dataType='drugScreen'").asDataFrame()


###pull files

singledose = []
multidose = []
for index,row in filelist.iterrows():
  #print(row['id'])
  dfile = pd.read_csv(syn.get(row['id']).path)
  dfile['improve_sample_id'] = row['specimenID']

  dfile = dfile.reset_index()
    
  #get single micromolar values
  sings = dfile[dfile.Concentration_uM==1]
  
  ##get counts of drugs
  dcounts = dfile.groupby("Drug").count().reset_index()
  more = dcounts[dcounts.Concentration_uM>1]['Drug']
  

  
  singledose.append(sings)
  multidose.append(dfile[dfile.Drug.isin(set(more))])


####first fit multidose curves...

fulltab = pd.concat(multidose)
#print(fulltab)
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
fulltab.to_csv('drug_response.tsv',sep='\t')


##fit curve
script='https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/build/utils/fit_curve.py'
subprocess.run(['wget',script])
subprocess.run(['python','fit_curve.py','--input','drug_response.tsv','--output','cnfDrugOutput'])


#####now we can take single drug points and format those
stab = pd.concat(singledose)

stab['time']=120
stab['time_unit']='hours'
stab['study']='cnfPDO'
stab['source']='synapse'
stab['improve_drug_id']=stab.Drug
stab['dose_response_value'] = stab.Viability_percentage/100.00
stab['dose_response_metric'] = 'uM_viability'

curve_cols = ['source','improve_sample_id','improve_drug_id','study','time','time_unit',
              'dose_response_metric','dose_response_value']


otab = pd.read_csv('cnfDrugOutput.0',sep='\t')
stab = stab[curve_cols]


newtab =pd.concat([otab,stab])

##rename file
newtab.to_csv('cohort1_curves.tsv',sep='\t',index=False)

##now add in single-point drug measurements

#store on synapse
syn.store(sc.File('cohort1_curves.tsv',parentId='syn65471813'))

