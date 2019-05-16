import pymongo, os, shutil
from pprint import pprint
from consensus import rgz_path, data_path, db, version
subjects = db['radio_subjects']

outdir = '/home/garon/Documents/RGZdata/sample'
os.mkdir(outdir)

sample = ['FIRSTJ102028.8+615133', 'FIRSTJ134340.3+635005', 'FIRSTJ174107.7+622523', 'FIRSTJ102104.7+641925', 'FIRSTJ134353.2+611451', 'FIRSTJ174107.9+622511', 'FIRSTJ102121.8+615540', 'FIRSTJ134354.4+634332', 'FIRSTJ174108.7+635605', 'FIRSTJ102123.6+632322', 'FIRSTJ134400.3+641857', 'FIRSTJ174113.8+630326', 'FIRSTJ102129.3+612711', 'FIRSTJ134402.8+620630', 'FIRSTJ174116.4+623742', 'FIRSTJ102130.3+611328', 'FIRSTJ134427.2+622814', 'FIRSTJ174122.5+633717', 'FIRSTJ102138.2+642613', 'FIRSTJ134434.3+641415']

for source in sample:
	subject = subjects.find_one({'metadata.source':source})
	zdir = outdir+'/'+subject['zooniverse_id']+'/'
	os.mkdir(zdir)
	shutil.copy(data_path+'/rgz/raw_images/RGZ-full.1/FIRST-IMGS/'+source+'.fits', zdir)
	for png in ['_radio.png', '_infrared.png', '_heatmap+contours.png']:
		shutil.copy(data_path+'/rgz/raw_images/RGZ-full.1/PNG-IMGS/'+source+png, zdir)
	shutil.copy(data_path+'/rgz/contours/'+str(subject['_id'])+'.json', zdir)

query = {'metadata.source':{'$in':sample}}
dump = 'mongodump -d radio -c radio_subjects -q "'+str(query)+'" -o '+outdir
os.system(dump.replace('$','\$'))