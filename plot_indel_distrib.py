# plot indel size distribution from mutatrix as boxplot

from fusion import read_vcf


df = read_vcf("tritici_60_genomes/mutatrix_60_genomes.vcf")
li = df["INFO"].str.split(";")
print(li[0])

li2 = [i[1] for i in li]
li3 = [i.split("=")[-1] for i in li2]
li4 = [i.split(",") for i in li3]
li4 = list(chain.from_iterable(li4))

print(li4[0])
print(len(li4))
# AC=3;LEN=38;NA=1;NS=60;TYPE=ins
li5 = [eval(i) for i in li4]


import plotly.express as px



## create some random data that looks similar
df = pd.DataFrame({'indel_size':li5})
fig = px.box(df, x="indel_size")

## loop through the values you want to label and add them as annotations
for x in zip(["min","q1","med","q3","max"],df.quantile([0,0.25,0.5,0.75,1]).iloc[:,0].values):
    fig.add_annotation(
        x=x[1],
        y=0.3,
        text=x[0] + ":" + str(x[1]),
        showarrow=False
        )
	

# fig = px.box(li5)
# fig.show()
fig.write_html('example_fig.html') 
# fig.write_image("fig1.png")