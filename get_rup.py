from libcomcat.search import get_event_by_id
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt

'''
Get rupture files of only unique earthquakes
other files containing this info?
'''

write_files = True 
use_big_eids = True 
path_to_rupture_files = '/Users/kksmith/im_scripts/data/rup/'

# # old one
# fname = '~/im_scripts/data/GMMout/NumRecs.txt'
# mydata = pd.read_csv(fname, sep='\t', header=None)
# UEIDs = mydata[0]

# proper one
fname = "~/Documents/GMM_data/for_kyle_conus/CONUS_HW_AK_h1_records.csv"
mydata = pd.read_csv(fname)

UEIDs, md_ue_loc = np.unique(mydata['EarthquakeId'], return_index=True)
ueq_data = mydata.loc[md_ue_loc,:]
#breakpoint()

# checking large EQs
bigeqs = ueq_data[ueq_data['EarthquakeMagnitude'] >= 6.5]
bigeqs['EarthquakeId']
Big_EIDs, be_loc = np.unique(bigeqs['EarthquakeId'], return_index=True)
event_pages = []
bad_connects = []
for b_i, bevid in enumerate(Big_EIDs):
    try:
        byid = get_event_by_id(eventid=bevid)
        event_pages.append(byid.toDict()['url'])
        preferred_product = byid.getProducts('shakemap')
        # if write_files:
        #     preferred_product[0].getContent('rupture.json','./data/rup/{}_rupture.json'.format(evid))
    except:
        bad_connects.append(bevid)
        print("connection issue w {} w/ M{}".format(bevid, bigeqs['EarthquakeMagnitude'][be_loc[b_i]]))
        continue

BC = pd.DataFrame(bad_connects)
BC.to_csv('bad_big_event_pages.csv', index=False, header=False)
EP = pd.DataFrame(event_pages)
EP.to_csv('big_event_pages.csv', index=False, header=False)
breakpoint()

if use_big_eids:
    EID_set = Big_EIDs 
    rup_tag = "_bigEQ_"
else:
    EID_set = UEIDs
    rup_tag = "_"

toteids = len(EID_set)
eids_w_rup = []
mags = []
eids_wO_rup = []
mags_no_rup = []
issue_eids = []
noeid = 0
ctr = 0

for uevid in EID_set:
    ctr += 1
    try:
        byid = get_event_by_id(eventid=uevid)
    except:
        issue_eids.append(uevid)
        continue
        # time.sleep(10)  
        # print("{}/{}".format(ctr,toteids))
        # byid = get_event_by_id(eventid=evid)

    try:
        preferred_product = byid.getProducts('shakemap')
        if write_files:
            preferred_product[0].getContent('rupture.json', path_to_rupture_files + '{}_rupture.json'.format(uevid))

        eids_w_rup.append(uevid)
        mags.append(byid['mag'])
        print('YES! Making rupture file of {} w/ M{}'.format(uevid,byid['mag']))
    except:
        noeid += 1
        print('NO! No shakemap product found for {} w/ M{}'.format(uevid,byid['mag']))

        eids_wO_rup.append(uevid)
        mags_no_rup.append(byid['mag'])

rup_eqs = pd.DataFrame(np.transpose([eids_w_rup,mags]))

if write_files:
    rup_eqs.to_csv("./data/avail" + rup_tag + "rup.txt", sep='\t', index=False, header=False)

no_rup_eqs = pd.DataFrame(np.transpose([eids_wO_rup,mags_no_rup]))
if write_files:
    no_rup_eqs.to_csv("./data/no" + rup_tag + "rup.txt", sep='\t', index=False, header=False)

print('missing {}/{} rupture files'.format(noeid,toteids))
print('{}/{} bad grabs'.format(len(issue_eids),toteids))
breakpoint()
print("end")
