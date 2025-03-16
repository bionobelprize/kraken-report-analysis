from process_dict import *
import sys,json,copy
import pandas as pd
from skbio.diversity import alpha_diversity
import plotkit as pk
def merge_dict(source,target):
	for key,item in source.items():
		if key in target:
			next_target=target[key]
			merge_dict(item,next_target)
		else:
			target.update({key:item})
			return
para_list=[(1,'samples'),(1,'subclade'),(3,'subclade'),(5,'subclade'),(7,'subclade'),(9,'subclade'),(11,'subclade'),(13,'subclade'),(15,'subclade'),(17,'subclade'),(19,'subclade'),(21,'subclade'),(23,'subclade')]
skip_list=[(3,'chromosomes')]
para_x_dict={}
target_c=('slevel','G')
locate_key=True
level=0
total_count=0
path='/home/zwb/yhd_sample.txt.combine.tsv.modified'
#path='/home/zwb/seawage_meta.sample.combine.tsv.modified'
tfile=open(path,'r')
ttdict=json.loads(tfile.read())
tfile.close()
aaa=dig_dict(level,ttdict,'',para_list,target_c,para_x_dict,total_count,locate_key,skip_list)
#print(len(aaa),type(aaa),type(aaa[0]),type(aaa[1]),type(aaa[2]),type(aaa[3]))
df=dict_to_df(aaa[1])
#filter virus reads
x1=pk.filter_by_c_v([df,'samples','target_species_name','target_percent',path],'7_species_name','cellular organisms')
x2=pk.normalized_percent(x1)
df=x2[0]
df=df[df['target_percent_norm']>0]
#al=pk.merge_parallel([df,'3_date','target_species_name','target_percent_norm'])
#print(al[1].index)
al_ser=df.set_index(['samples','target_species_name'])['target_percent_norm']
data,group_l=pk.matrix_prepare(al_ser)
adiv_sobs=alpha_diversity('sobs',data,group_l)
adiv_sobs_dict=adiv_sobs.to_dict()
recover_dict={'samples':{}}
for x,y in adiv_sobs_dict.items():
	recover_dict['samples'].update({x:{'alpha_diversity':y}})
merge_dict(recover_dict,ttdict)
js = json.dumps(ttdict,indent=4)
efile = open(path+'.test','w')
efile.write(js)
efile.close()