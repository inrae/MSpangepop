import pandas as pd
import yaml, random, re

## read yaml file containing variants informations
def read_yaml(f):
	stream = open(f, 'r')
	d = yaml.safe_load(stream)
	stream.close()
	if sum(d.values()) != 100:
		raise ValueError("Sum of values must be 100")
	return(d)

infos_dict = {'deletion': None, 'insertion': None, 'inversion': None, 'SNP': None,
		 'tandem duplication': 2, 'inverted tandem duplication': 2, 
		 'translocation copy-paste': 0, 'translocation cut-paste': 0, 'reciprocal translocation': 0}

## generate info field for translocations
# randomly select chromosome to add sequence to
def select_chr(fai):
	df = pd.read_table(fai, header=None, usecols=[0,1], names =["name", "length"])
	names=df["name"].to_list()
	chr = random.choice(names)
	df = pd.pivot_table(df, columns='name')
	len = df.iloc[0][chr]
	pos = random.randrange(len)
	li = [chr, pos]
	return(li)

# randomly select sequence orientation
def select_orientation(li):
	li.append(random.choice(["forward", "reverse"]))
	return(li)

## create dataframe with all the variants to simulate
def generate_type(nvar, yml, fai):
	d = read_yaml(yml)
	print(d)
	## list for variant type and info
	vartype = []
	infos = []
	## filter dict to keep value != 0
	dict = {k: v for k, v in d.items() if v != 0}
	for k in dict:
		# generate n variants
		perc = dict[k]
		n = round(perc*nvar/100)
		# add info field
		if re.search("translocation", k):
			for i in range(n):
				info = select_chr(fai)
				info = select_orientation(info)		
				if k == "reciprocal translocation":
					info = select_orientation(info)
				infos.extend([info])
				i += 1
		else:
			info = infos_dict[k]
			infos.extend([info]*n)
		vartype.extend([k]*n)
		
	for_df = {'variant': vartype, 'info': infos}

	df = pd.DataFrame(data = for_df)
	# shuffle values in dataframe
	df = df.sample(frac = 1)
	return(df)