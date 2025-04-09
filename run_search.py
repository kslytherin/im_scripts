"just for playing around with libcomcat"

from datetime import datetime
from libcomcat.search import search
from libcomcat.search import get_event_by_id
from libcomcat.dataframes import get_detail_data_frame
from libcomcat.dataframes import get_summary_data_frame
from libcomcat.dataframes import _describe_finite_fault

path_to_rupture_files = '/Users/kksmith/data/rup/'

bysearch = search(starttime=datetime(1994, 1, 17, 12, 30), endtime=datetime(1994, 4, 18, 12, 35), updatedafter=datetime(1994, 2, 18, 12, 35), producttype='shakemap', minmagnitude=4, maxmagnitude=8)
bys_frame = get_detail_data_frame(bysearch)
get_summary_data_frame(bysearch)

evid = 'ci3144585'

byid = get_event_by_id(eventid=evid, includedeleted=True)
preferred_product = byid.getProducts('shakemap')
preferred_product[0].contents
preferred_product[0].getContent('rupture.json', path_to_rupture_files + '{}_rupture.json'.format(evid))

breakpoint()


#product_type = search(starttime=datetime(1994, 1, 17, 12, 30), endtime=datetime(1994, 4, 18, 12, 35),producttype='shakemap')
breakpoint()
print("end")
