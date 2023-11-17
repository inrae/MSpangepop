import pandas as pd
import yaml, random

def read_yaml(f):
	stream = open(f, 'r')
	d = yaml.safe_load(stream)
	stream.close()
	if sum(d.values()) != 100:
		raise ValueError("Sum of values must be 100")

	return(d)

infos = {'deletion': None, 'insertion': None, 'inversion': None, 
		 'tandem duplication': 2, 'inverted tandem duplication': 2, 
		 'translocation copy-paste': 0, 'translocation cut-paste': 0, 'reciprocal translocation': 0, 'SNP': None}

# generate info field
def select_chr(fai):
	df = pd.read_table(fai, header=None, usecols=[0,1], names =["name", "length"])
	names=df["name"].to_list()
	chr = random.choice(names)
	df = pd.pivot_table(df, columns='name')
	len = df.iloc[0][chr]
	pos = random.randrange(len)
	li = [chr, pos]
	return(li)

def select_orientation(li):
	li.append(random.choice(["forward", "reverse"]))
	return(li)


def generate_type(nvar, yml, fai):
	d = read_yaml(yml)
	print(d)
	## list for variant type and info
	vartype = []
	infos = []
	## filter dict to keep value != 0
	dict = {k: v for k, v in d.items() if v != 0}
	for k in dict:
		perc = dict[k]
		n = round(perc*nvar/100)
		if k == "tandem duplication":
			info = 2
		elif k == 'translocation copy-paste' or k == 'translocation cut-paste':
			info = select_chr(fai)
			info = select_orientation(info)
		elif k == "reciprocal translocation":
			info = select_chr(fai)
			info = select_orientation(info)
			info = select_orientation(info)
		else:
			info = None

		vartype.extend([k]*n)
		infos.extend([info]*n)
	for_df = {'variant': vartype, 'info': infos}

	df = pd.DataFrame(data = for_df)
	df = df.sample(frac = 1)
	print(df)


# tandem duplication: 0
# inverted tandem duplication: 0
# translocation copy-paste: 0
# translocation cut-paste : 0
# reciprocal translocation: 0

	

fai = "/home/sukanya/tests/02_data/Sibirica_v1.0.fa.fai"
# select_chr(fai)
generate_type(12, "visor_sv_type.yaml", fai)
